function add_base_pair(tag1, tag2, edge1, edge2, LW_orientation, mode)
% add_base_pair(tag1, tag2, edge1, edge2, LW_orientation, mode)
%
% 在当前绘图 (gca) 上交互式地添加或删除碱基对 linker。
%
% 使用示例：
%   % 1) 添加一个非典型配对（或茎内配对，视当前 helix/stem 信息而定）
%   %    A:A:10 与 A:A:45 之间的 W/W/C 配对：
%   add_base_pair('A:A:10','A:A:45','W','W','C');
%
%   % 2) 删除这两个残基之间的所有 base-pair 类型 linker：
%   add_base_pair('A:A:10','A:A:45',[],[],[],'delete');
%
% 说明：
%   - tag1, tag2 使用与 FASTA / .base_pairs.txt 一致的形式：
%       'A:10'、'A:A:10'、'B:QA:25' 等（可带或不带 segid）。
%   - edge1, edge2 分别是第一个、第二个残基的边缘类型：'W'/'H'/'S'。
%   - LW_orientation 为配对几何类型：'C' 或 'T'。
%   - mode:
%       * 省略或为 'add' 时：添加（或覆盖）碱基对 linker；
%       * 为 'delete'/'remove' 时：删除这两个残基之间的 base-pair linker。
%
% 注意：
%   - 本函数依赖于当前 axes (gca) 中已经由 initialize_drawing 等函数创建的
%     Residue_* 与 Helix_* 对象。
%   - 添加时内部会调用 setup_base_pair_linkers 和 redraw_helices，以确保新的
%     linker 立即可见。
%
% (C) R. Das, Stanford University, 2017
% 增强与中文注释 (C) OpenAI 协助整理, 2025

if ~exist('mode','var') || isempty(mode)
    mode = 'add';
end

if nargin < 3 || isempty(edge1);          edge1 = 'W'; end
if nargin < 4 || isempty(edge2);          edge2 = 'W'; end
if nargin < 5 || isempty(LW_orientation); LW_orientation = 'C'; end

mode = lower(mode);

switch mode
    case {'add','create'}
        add_base_pair_internal(tag1, tag2, edge1, edge2, LW_orientation);
    case {'delete','remove'}
        delete_base_pair_internal(tag1, tag2);
    otherwise
        error('add_base_pair: 未识别的 mode = %s，应为 ''add'' / ''delete''。', mode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function add_base_pair_internal(tag1, tag2, edge1, edge2, LW_orientation)

fprintf('[add_base_pair] 尝试添加碱基对：%s - %s，edges=%s/%s，LW=%s\n', ...
    tag1, tag2, edge1, edge2, LW_orientation);

% 解析两个残基标签
[resnum1, chain1, segid1, ok1] = get_one_resnum_from_tag(tag1);
[resnum2, chain2, segid2, ok2] = get_one_resnum_from_tag(tag2);

if ~ok1
    error('[add_base_pair] 无法解析第一个残基标签：%s', tag1);
end
if ~ok2
    error('[add_base_pair] 无法解析第二个残基标签：%s', tag2);
end

% 构造 base_pair 结构体
clear base_pair
base_pair.resnum1 = resnum1;
base_pair.chain1  = chain1;
base_pair.segid1  = segid1;
base_pair.resnum2 = resnum2;
base_pair.chain2  = chain2;
base_pair.segid2  = segid2;
base_pair.edge1   = edge1;
base_pair.edge2   = edge2;
base_pair.LW_orientation = upper(LW_orientation);

% 规范化顺序（按 chain/segid/resnum 排序）
base_pair = ordered_base_pair(base_pair);

% 在当前绘图中安装该 base_pair 对应的 linker
setup_base_pair_linkers({base_pair});

% 重新绘制所有 helices，以确保新 linker 被绘制出来
fprintf('[add_base_pair] 已安装 linker，重新绘制所有 helices 以显示新的碱基对。\n');
redraw_helices();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delete_base_pair_internal(tag1, tag2)

fprintf('[add_base_pair] 尝试删除碱基对 linker：%s - %s\n', tag1, tag2);

% 解析两个残基标签
[resnum1, chain1, segid1, ok1] = get_one_resnum_from_tag(tag1);
[resnum2, chain2, segid2, ok2] = get_one_resnum_from_tag(tag2);

if ~ok1
    error('[add_base_pair] 无法解析第一个残基标签：%s', tag1);
end
if ~ok2
    error('[add_base_pair] 无法解析第二个残基标签：%s', tag2);
end

res_tag1 = sanitize_tag(sprintf('Residue_%s%s%d', chain1, segid1, resnum1));
res_tag2 = sanitize_tag(sprintf('Residue_%s%s%d', chain2, segid2, resnum2));

if ~isappdata(gca, res_tag1)
    fprintf('[add_base_pair] 警告：当前图中找不到残基 %s，对应 tag=%s。\n', tag1, res_tag1);
    return;
end
if ~isappdata(gca, res_tag2)
    fprintf('[add_base_pair] 警告：当前图中找不到残基 %s，对应 tag=%s。\n', tag2, res_tag2);
    return;
end

residue1 = getappdata(gca, res_tag1);

if ~isfield(residue1, 'linkers') || isempty(residue1.linkers)
    fprintf('[add_base_pair] 残基 %s 当前没有任何 linker，无法删除 base-pair linker。\n', res_tag1);
    return;
end

% 要删除的 linker 类型（base pair 相关）
basepair_types = {'stem_pair','noncanonical_pair','long_range_stem_pair'};

linkers_to_delete = {};

for i = 1:length(residue1.linkers)
    linker_tag = residue1.linkers{i};
    if ~isappdata(gca, linker_tag)
        continue;
    end
    linker = getappdata(gca, linker_tag);
    if ~isfield(linker, 'residue1') || ~isfield(linker, 'residue2')
        continue;
    end

    % 检查是否恰好连在这两个残基之间（顺序不敏感）
    pair_match = (strcmp(linker.residue1, res_tag1) && strcmp(linker.residue2, res_tag2)) || ...
                 (strcmp(linker.residue1, res_tag2) && strcmp(linker.residue2, res_tag1));
    if ~pair_match
        continue;
    end

    % 仅删除 base-pair 类型的 linker
    if ~isfield(linker, 'type') || ~any(strcmp(linker.type, basepair_types))
        fprintf('[add_base_pair] 跳过非 base-pair 类型的 linker：%s (type=%s)\n', ...
            linker_tag, getfield(linker, 'type', 'UNKNOWN'));
        continue;
    end

    fprintf('[add_base_pair] 标记待删除 linker：%s (type=%s)\n', linker_tag, linker.type);
    linkers_to_delete = [linkers_to_delete, {linker_tag}];
end

if isempty(linkers_to_delete)
    fprintf('[add_base_pair] 未找到需要删除的 base-pair linker。\n');
    return;
end

% 删除所有匹配的 linker
delete_linker(linkers_to_delete, 1);
fprintf('[add_base_pair] 已删除 %d 个 base-pair linker。\n', length(linkers_to_delete));

