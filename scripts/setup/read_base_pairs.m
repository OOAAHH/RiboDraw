function base_pairs = read_base_pairs( base_pairs_file, resnum_ref, chains_ref, segid_ref )
% base_pairs = read_base_pairs( base_pairs_file, resnum_ref, chains_ref, segid_ref )
%
% 读取由Rosetta rna_motif程序输出的.base_pairs.txt文件
% 该文件包含所有Watson-Crick和非Watson-Crick碱基对信息
%
% 输入参数:
%  base_pairs_file = 文本文件，每行格式如下:
%                    A:1  B:20 W H C
%                    即:
%                    chain1[:segid1]:resnum1 chain2[:segid2]:resnum2  edge1  edge2  LW_orientation
%
%  可选参数（用于更鲁棒地补全 segid）:
%   resnum_ref, chains_ref, segid_ref  = 来自 FASTA 的索引信息，长度应相同，
%                                        分别为残基编号、链 ID、段 ID。只有在这三个
%                                        参数都提供且非空时，才会尝试为缺失 segid 的
%                                        basepair 自动补全段 ID。
%
%                    其中:
%                    edge1, edge2 = W/H/S (Watson-Crick/Hoogsteen/Sugar边缘)
%                    LW_orientation = C/T (顺式/反式，基于两个碱基的糖苷键相对于碱基间连接器的方向)
%
% 输出参数:
%  base_pairs = 结构体数组，包含相同的信息，但经过重新排序
%               排序规则：链/段ID较小的残基在前，如果相同则残基编号较小的在前
%
% 相关函数: SETUP_BASE_PAIRS, READ_BASE_STACKS
% 
% (C) R. Das, Stanford University, 2017

% 初始化空数组来存储碱基对信息
base_pairs = {};

% 显示当前工作目录和文件路径，用于调试
disp(['[read_base_pairs] 当前工作目录: ', pwd]);
disp(['[read_base_pairs] 尝试打开文件: ', base_pairs_file]);

% 判断是否启用基于 FASTA 信息的 segid 补全逻辑
use_segid_fix = false;
segid_map = [];
fixed_segid_count = 0;
missing_segid_count = 0;

if exist('resnum_ref','var') && exist('chains_ref','var') && exist('segid_ref','var') && ...
        ~isempty(resnum_ref) && ~isempty(chains_ref) && ~isempty(segid_ref)
    use_segid_fix = true;
    fprintf('[read_base_pairs] 收到参考序列信息，将尝试为缺失 segid 的碱基对补全段 ID。\n');
    assert(length(resnum_ref) == length(chains_ref), ...
        '[read_base_pairs] resnum_ref 与 chains_ref 长度不一致');
    assert(length(resnum_ref) == length(segid_ref), ...
        '[read_base_pairs] resnum_ref 与 segid_ref 长度不一致');

    % 构建 (chain,resnum) -> segid 的映射，只记录非空 segid
    segid_map = containers.Map();
    for i = 1:length(resnum_ref)
        key = sprintf('%s:%d', chains_ref(i), resnum_ref(i));
        segid_val = segid_ref{i};
        if isempty(segid_val)
            continue;
        end
        if isKey(segid_map, key)
            if ~strcmp(segid_map(key), segid_val)
                fprintf(['[read_base_pairs] 警告：参考序列中 key=%s 存在多个 segid 候选：' ...
                         '''%s'' 与 ''%s''，保留第一个。\n'], ...
                        key, segid_map(key), segid_val);
            end
        else
            segid_map(key) = segid_val;
        end
    end
    fprintf('[read_base_pairs] segid 映射构建完成，非空 segid 条目数 = %d。\n', segid_map.Count);
else
    fprintf('[read_base_pairs] 未提供参考序列信息，按原逻辑读取 base_pairs（不尝试 segid 修复）。\n');
end

% 确保base_pairs_file是字符向量
if isstring(base_pairs_file)
    base_pairs_file = char(base_pairs_file);
end

% 检查文件是否存在，如果不存在则直接返回
if ~exist( base_pairs_file, 'file' )
    disp(['错误：文件不存在: ', base_pairs_file]);
    return;
end

% 打开文件
fid = fopen( base_pairs_file );
if fid == -1
    disp(['[read_base_pairs] 错误：无法打开文件: ', base_pairs_file]);
    return;
end

% 逐行读取文件内容
while ~feof( fid )
    % 读取一行
    line = fgetl( fid );
    % 按空格分割行内容
    cols = strsplit( line, ' ' );
    
    % 确保至少有5列数据
    if length( cols ) >= 5        
        % 从第一列解析第一个残基的信息（残基号、链ID、段ID）
        [base_pair.resnum1,base_pair.chain1,base_pair.segid1] = get_one_resnum_from_tag( cols{1} );
        % 从第二列解析第二个残基的信息
        [base_pair.resnum2,base_pair.chain2,base_pair.segid2] = get_one_resnum_from_tag( cols{2} );

        % 如果启用了 segid 补全逻辑，则尝试为缺失 segid 的残基补全段 ID
        if use_segid_fix
            for which = 1:2
                segid_field = sprintf('segid%d', which);
                chain_field = sprintf('chain%d', which);
                res_field   = sprintf('resnum%d', which);

                if isempty(base_pair.(segid_field))
                    key = sprintf('%s:%d', base_pair.(chain_field), base_pair.(res_field));
                    if isKey(segid_map, key)
                        new_segid = segid_map(key);
                        fprintf(['[read_base_pairs] 为 %s 修复 segid：链 %s 残基 %d ' ...
                                 '-> segid ''%s''（来自 FASTA 映射）。\n'], ...
                                cols{which}, base_pair.(chain_field), base_pair.(res_field), new_segid);
                        base_pair.(segid_field) = new_segid;
                        fixed_segid_count = fixed_segid_count + 1;
                    else
                        % 找不到匹配项，记录一次缺失，方便之后汇总
                        missing_segid_count = missing_segid_count + 1;
                        fprintf(['[read_base_pairs] 警告：无法为 %s（链 %s 残基 %d）找到 segid，' ...
                                 '对应 FASTA 中无匹配项。\n'], ...
                                cols{which}, base_pair.(chain_field), base_pair.(res_field));
                    end
                end
            end
        end

        % 存储边缘信息
        base_pair.edge1 = cols{3};
        base_pair.edge2 = cols{4};
        % 存储LW方向信息
        base_pair.LW_orientation = cols{5};
        % 将解析后的碱基对信息添加到结果数组中，并确保正确的排序
        base_pairs = [base_pairs,ordered_base_pair(base_pair)];
    end;
end

% 关闭文件
fclose(fid);

% 打印 segid 修复统计信息
if use_segid_fix
    fprintf('[read_base_pairs] segid 修复完成：成功补全 %d 个残基的 segid，仍有 %d 个残基 segid 为空。\n', ...
        fixed_segid_count, missing_segid_count);
end

