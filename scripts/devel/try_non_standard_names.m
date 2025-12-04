function try_non_standard_names( sequence, resnum, chains, segid, non_standard_residues );
% try_non_standard_names( sequence, resnum, chains, non_standard_residues);
%
% 此函数用于处理非标准残基的命名，
% 将非标准残基的名称转换为首选显示名称，并更新对应的 residue 对象。
%
% 输入参数：
%  sequence                - 序列字符串
%  resnum                  - 每个序列位置对应的残基编号 (整数)
%  chains                  - 每个序列位置对应的链标识 (字符)
%  segid                   - 每个序列位置对应的段标识 (字符串，可能为空)
%  non_standard_residues   - 结构体，包含非标准残基的位置索引 (index) 与其对应的名称 (name)
%
% (C) R. Das, Stanford University

for i = 1:length(sequence)
    % 提取当前残基的链、残基编号和段标识
    chain = chains(i);
    res   = resnum(i);
    seg   = segid{i};
    
    % 构造当前残基的标签，格式为 'Residue_{chain}{seg}{res}'
    res_tag = sanitize_tag(sprintf('Residue_%s%s%d', chain, seg, res));
    fprintf('【调试】处理第 %d 个残基：chain = %s, res = %d, seg = %s, res_tag = %s\n', i, chain, res, seg, res_tag);
    
    % 检查当前坐标轴 (gca) 中是否存在对应的应用数据，如果不存在则跳过该残基
    if ~isappdata(gca, res_tag)
        fprintf('【警告】当前坐标轴中不存在标签 %s 的应用数据，跳过该残基。\n', res_tag);
        continue;
    end;
    
    % 获取当前残基对应的 residue 对象
    residue = getappdata(gca, res_tag);
    fprintf('【调试】当前的 residue 对象字段：\n');
    disp(fieldnames(residue));
    if isfield(residue, 'helix_tag')
        fprintf('【调试】当前 residue 的 helix_tag: %s\n', residue.helix_tag);
    else
        fprintf('【调试】当前 residue 不包含 helix_tag 字段。\n');
    end

    
    % 计算当前残基在序列中的位置 (seqpos)
    % 条件：链匹配、残基编号相等且 segid 匹配
    seqpos = intersect(intersect(strfind(chains, chain), find(resnum == res)), ...
                       find(strcmp(segid, seg)));
    fprintf('【调试】计算得到 seqpos: %s\n', mat2str(seqpos));
    
    % 更新 residue 的名称信息：使用对应位置的 sequence 字符转换为大写作为显示名称
    residue.name = upper(sequence(seqpos));
    residue.original_name = sequence(seqpos);
    fprintf('【调试】更新后 residue.name: %s, original_name: %s\n', residue.name, residue.original_name);
    
    % 检查当前序列位置是否为非标准残基
    if any(non_standard_residues.index == seqpos)
        idx = find(non_standard_residues.index == seqpos);
        name = non_standard_residues.name{idx};
        fprintf('【调试】发现非标准残基：位置 %s 对应名称: %s\n', mat2str(seqpos), name);
        
        % 保存非标准残基名称，并更新 residue.name 为首选显示名称
        residue.non_standard_residue_name = name;
        residue.name = get_preferred_display_name(name);
        fprintf('【调试】更新后的 residue.name (非标准): %s\n', residue.name);
        
        % 如果更新后的名称为 'X'，则输出警告信息
        if strcmp(residue.name, 'X')
            fprintf('【警告】残基显示名称为 X，原始非标准名称为 %s\n', name);
        end;
    end
    
    % 如果 residue 对象包含 'handle' 字段，则更新对应图形对象的显示字符串
    if isfield(residue, 'handle')
        set(residue.handle, 'String', residue.name);
        fprintf('【调试】更新图形对象显示为: %s\n', residue.name);
    end;
    
    % 将更新后的 residue 对象存回当前坐标轴 (gca) 的应用数据中
    setappdata(gca, res_tag, residue);
    fprintf('【调试】存储更新后的 residue 对象到 gca，标签为: %s\n', res_tag);
    fprintf('------------------------------------------------------\n');
end
