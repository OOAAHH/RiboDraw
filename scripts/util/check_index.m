function check_index(base_pairs, resnum, chains)
% check_index 检查 base_pairs 与 FASTA 文件中残基索引的一致性
% check_index verifies that the residue indices extracted from 
%
% base_pairs match the indices generated from the FASTA file.
%
% 用法/Usage:
%   check_index(base_pairs, resnum, chains)
%
% 输入/Inputs:
%   base_pairs - cell 数组，每个单元为一个结构体，包含字段：
%                resnum1, chain1, segid1, resnum2, chain2, segid2, ... 等。
%
%                English: a cell array where each element is a struct with fields 
%                such as resnum1, chain1, segid1, resnum2, chain2, segid2, etc.
%
%   resnum     - 数组，包含从 FASTA 文件中读取的所有残基编号。
%
%                English: a numeric array containing residue numbers extracted 
%                from the FASTA file.
%
%   chains     - 字符串或字符数组，每个元素对应 FASTA 文件中同一位置的链标识。
%
%                English: a string or character array where each element 
%                corresponds to the chain identifier for a residue in the FASTA file.
%
% 功能描述/Description:
%   从 base_pairs 中提取残基索引（例如 "A,1"、"B,16"），并将其与由 FASTA 文件中的 
%   resnum 和 chains 生成的索引进行比较。如果存在不匹配，则输出警告信息。
%
%                English: The function extracts residue indices (e.g., "A,1", "B,16")
%                from base_pairs and compares them with the indices generated from
%                the FASTA file using resnum and chains. A warning is issued if any
%                mismatch is found.
%
% (C) Hao Sun, GuangZhou National Labortory, 2025

%% 1. 从 base_pairs 中提取残基索引 / Extract indices from base_pairs
bp_indices = {};  % 用于存储 base_pairs 中提取的 "chain,resnum" 字符串
% English: bp_indices stores the "chain,resnum" strings extracted from base_pairs.
fprintf('【调试/Debug】开始遍历 base_pairs...\n');

for i = 1:length(base_pairs)
    bp = base_pairs{i};
    % 提取第一个残基信息：chain1 和 resnum1
    % English: Extract first residue info: chain1 and resnum1.
    idx1 = sprintf('%s,%d', bp.chain1, bp.resnum1);
    bp_indices{end+1} = idx1;
    fprintf('【调试/Debug】 base_pairs{%d} -> idx1: %s\n', i, idx1);
    
    % 提取第二个残基信息：chain2 和 resnum2
    % English: Extract second residue info: chain2 and resnum2.
    idx2 = sprintf('%s,%d', bp.chain2, bp.resnum2);
    bp_indices{end+1} = idx2;
    fprintf('【调试/Debug】 base_pairs{%d} -> idx2: %s\n', i, idx2);
end

%% 2. 生成 FASTA 文件中的残基索引 / Generate FASTA indices from resnum and chains
fasta_indices = cell(1, length(resnum));
for i = 1:length(resnum)
    % 假设 chains 中的每个字符对应一个残基
    % English: Assume each element in chains corresponds to one residue.
    fasta_indices{i} = sprintf('%s,%d', chains(i), resnum(i));
end

fprintf('【调试/Debug】 FASTA 索引生成完成，共 %d 个残基。\n', length(fasta_indices));
% 如有需要，可打印所有 FASTA 索引
% disp(fasta_indices);

%% 3. 比较 base_pairs 索引与 FASTA 索引 / Compare base_pairs indices with FASTA indices
inconsistency_found = false;
for i = 1:length(bp_indices)
    idx = bp_indices{i};
    if ~any(strcmp(fasta_indices, idx))
        warning('【警告/Warning】base_pairs 中的索引 %s 在 FASTA 中未找到！', idx);
        inconsistency_found = true;
    else
        fprintf('【调试/Debug】 索引 %s 在 FASTA 中存在。\n', idx);
    end
end

if ~inconsistency_found
    fprintf('【调试/Debug】 所有 base_pairs 索引均在 FASTA 中找到，一致性检查通过。\n');
end
