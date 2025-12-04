function check_index(base_pairs, resnum, chains)
% check_index 检查 base_pairs 中的残基信息与 FASTA 文件中的残基索引是否一致
%
% 调用方式：
%    check_index(base_pairs, resnum, chains)
%
% 输入参数：
%   base_pairs - cell 数组，每个单元为一个结构体，包含字段：
%                resnum1, chain1, segid1, resnum2, chain2, segid2, ...等，
%                例如：
%                  base_pairs{1}.resnum1 = 1;
%                  base_pairs{1}.chain1  = 'A';
%                  base_pairs{1}.resnum2 = 16;
%                  base_pairs{1}.chain2  = 'B';
%
%   resnum   - 数组，包含从 FASTA 中读取的所有残基编号，
%              例如：[1 2 3 4 ... 15 1 2 3 ... 16]
%
%   chains   - 字符串或字符数组，每个元素对应 resnum 中同一位置的链标识，
%              例如：'AAAAAAAAAAAAAAABBBBBBBBBBBBBBBB'
%
% 该函数会将 base_pairs 中每个配对提取的残基（以"链,残基编号"形式，如 "A,1"）与
% 由 FASTA 中的 chains 和 resnum 生成的索引逐一比较。如果存在 base_pairs 中的某个
% "链,残基编号"未能在 FASTA 索引中找到，则输出警告提示不一致。
%
% (C) R. Das, Stanford University.

%% 1. 从 base_pairs 中提取残基索引
bp_indices = {};  % 存放 base_pairs 中提取的 "chain,resnum" 字符串

fprintf('【调试】开始遍历 base_pairs...\n');
for i = 1:length(base_pairs)
    bp = base_pairs{i};
    % 提取第一个残基信息：chain1 和 resnum1
    idx1 = sprintf('%s,%d', bp.chain1, bp.resnum1);
    bp_indices{end+1} = idx1;
    fprintf('【调试】 base_pairs{%d} -> idx1: %s\n', i, idx1);
    
    % 提取第二个残基信息：chain2 和 resnum2
    idx2 = sprintf('%s,%d', bp.chain2, bp.resnum2);
    bp_indices{end+1} = idx2;
    fprintf('【调试】 base_pairs{%d} -> idx2: %s\n', i, idx2);
end

%% 2. 生成 FASTA 中的残基索引
fasta_indices = cell(1, length(resnum));
for i = 1:length(resnum)
    % 这里假设 chains 是一个字符数组或字符串，其中每个字符对应一个残基
    fasta_indices{i} = sprintf('%s,%d', chains(i), resnum(i));
end

fprintf('【调试】 FASTA 索引生成完成，共 %d 个残基。\n', length(fasta_indices));
% 如有需要，可打印全部 FASTA 索引
% disp(fasta_indices);

%% 3. 检查 base_pairs 中的每个索引是否都存在于 FASTA 索引中
inconsistency_found = false;
for i = 1:length(bp_indices)
    idx = bp_indices{i};
    if ~any(strcmp(fasta_indices, idx))
        warning('basepairs的信息和fasta不一致：basepairs中的索引 %s 在fasta中未找到！', idx);
        inconsistency_found = true;
    else
        fprintf('【调试】 索引 %s 在 FASTA 中存在。\n', idx);
    end
end

if ~inconsistency_found
    fprintf('【调试】 所有 base_pairs 索引均在 FASTA 中找到，一致性检查通过。\n');
end
