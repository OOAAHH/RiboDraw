function initialize_drawing(tag )
% initialize_drawing( tag )
%
% 主函数，用于启动 RiboDraw 绘图工具，基于 Rosetta rna_motif 运行 PDB 文件的输出结果。
%
% 输入参数：
%   tag = 输入的 PDB 文件名称，可能包含路径
%         例如，如果 tag 为 "RNA.pdb"，则该函数期望存在以下文件：
%              RNA.pdb.fasta
%              RNA.pdb.base_pairs.txt
%              RNA.pdb.stacks.txt
%              RNA.pdb.other_contacts.txt
%              RNA.pdb.stems.txt
%              RNA.pdb.ligands.txt
%              RNA.pdb.motifs.txt
%
% (C) R. Das, Stanford University, 2017

% 添加调试信息
disp(['tag的类型: ', class(tag)]);
disp(['tag的内容: ', tag]);
disp(['tag的长度: ', num2str(length(tag))]);

%% 1. 读取输入文件数据
% -------------------------
% 从 .fasta 文件中获取序列、残基编号、链信息、段标识以及非标准残基信息
disp("1. 开始读取信息");

% 确保tag是字符向量，兼容所有MATLAB版本
if isstring(tag)
    % 如果是string类型（R2017a及以后版本），转换为char
    tag = char(tag);
end

% 准备文件名
fasta_file = [tag, '.fasta'];
base_pairs_file = [tag, '.base_pairs.txt'];
stacks_file = [tag, '.stacks.txt'];
other_contacts_file = [tag, '.other_contacts.txt'];
stems_file = [tag, '.stems.txt'];
ligands_file = [tag, '.ligands.txt'];
motifs_file = [tag, '.motifs.txt'];

% 显示文件名信息
disp(['fasta文件: ', fasta_file]);
disp(['base_pairs文件: ', base_pairs_file]);

[sequence, resnum, chains, segid, non_standard_residues] = get_sequence(fasta_file);
disp(tag);
% 读取碱基对信息，包含非典型配对
base_pairs = read_base_pairs(base_pairs_file, resnum, chains, segid);

% 检查索引信息和fasta中索引信息是否一致
check_index(base_pairs, resnum, chains)

% 读取碱基堆叠信息，包含非典型堆叠
base_stacks = read_base_stacks(stacks_file);

% 读取其它接触信息（可能包括非标准的分子间接触）
other_contacts = read_other_contacts(other_contacts_file);

% 读取 RNA 二级结构中茎（helix）部分的信息
stems = read_stems(stems_file);

% 读取配体信息，如小分子、离子等
ligands = read_ligands(ligands_file);

% 读取结构模体（motifs）信息
motifs = read_motifs(motifs_file);

% 如果碱基对信息为空但茎信息存在，则根据茎信息生成碱基对
if length(base_pairs) == 0 & length(stems) > 0
    disp("碱基对信息为空但茎信息存在，根据茎信息生成碱基对");
    base_pairs = get_base_pairs_from_stems(stems);
end
disp("1. 读取完成");

%% 2. 设置绘图环境
% -------------------------
% 清空当前图形窗口，并设置当前坐标轴（gca）的大小为全屏
disp("2. 设置绘图环境");
clf; 
set(gca, 'Position', [0 0 1 1]);

% 保持当前绘图状态，允许多次绘制
hold on;

% 初始化一个与序列长度相同的零向量（可能用于后续计算或标记）
t = zeros(1, length(sequence));

% 设置绘图区域的坐标范围
axis([0 200 0 200]);

% 将默认的绘图设置保存到当前坐标轴的应用数据中，便于后续各函数调用时使用
setappdata(gca, 'plot_settings', default_plot_settings);
disp("2. 绘图环境设置完成");

%% 3. 初始化结构对象及数据
% -------------------------
% 设置茎（helix）对象的默认位置，例如确定每个螺旋的中心位置
disp("3. 初始化结构对象及数据，开始处理helix");
disp("3. 设置茎（helix）对象的默认位置");
stems = set_default_stem_positions(stems);

% 根据序列、残基编号、链信息和段标识初始化每个残基对象，同时设置每个残基所属的 helix_tag
disp("3.  根据序列、残基编号、链信息和段标识初始化每个残基对象，同时设置每个残基所属的 helix_tag");
stems = setup_residues(stems, sequence, resnum, chains, segid);

% 建立茎之间的配对关系（例如，为后续绘制联系线做准备）
disp("3. 设置茎（helix）之间的配对关系");
setup_stem_partner(stems);

%% 4. 初始化连接器（Linkers）
% -------------------------
% 设置箭头连接器，通常用于表示方向性或序列连接，依据残基编号、链信息和段标识
disp("4. 初始化连接器");
disp("4. 设置箭头连接器，通常用于表示方向性或序列连接，依据残基编号、链信息和段标识");
setup_arrow_linkers(resnum, chains, segid);

% 根据碱基对信息设置连接器，用于在图中显示碱基对之间的联系
disp("4. 根据碱基对信息设置连接器，用于在图中显示碱基对之间的联系");
setup_base_pair_linkers(base_pairs);

% 根据碱基堆叠信息设置连接器，用于显示碱基堆叠的相互作用
disp("4. 根据碱基堆叠信息设置连接器，用于显示碱基堆叠的相互作用");
setup_base_stack_linkers(base_stacks);

% 根据其它接触信息设置连接器，用于显示非典型接触
disp("4. 根据其它接触信息设置连接器，用于显示非典型接触");
setup_other_contact_linkers(other_contacts);

% 尝试处理非标准残基的名称，可能进行名称规范化或转换
disp("4. 尝试处理非标准残基的名称，可能进行名称规范化或转换");
try_non_standard_names(sequence, resnum, chains, segid, non_standard_residues);

%% 5. 初始化其它结构（模体、配体等）
% -------------------------
% 以下两行代码（关于同轴堆叠）被注释掉，如果需要可以启用：
% coaxial_stacks = get_coaxial_stacks(base_pairs, base_stacks, stems);
% setup_coaxial_stacks(coaxial_stacks);
disp("5. 初始化其它结构（模体、配体等）");
% 设置配体信息，准备绘制配体分子
disp("5. 设置配体信息，准备绘制配体分子");
setup_ligands(ligands);

% 设置模体信息，准备绘制结构模体
disp("5. 设置模体信息，准备绘制结构模体");
setup_motifs(motifs);

%% 6. 绘制结构及后续设置
% -------------------------
% 绘制 RNA 螺旋（helix）结构，根据茎信息生成图形表示
disp("6. 绘制 RNA 螺旋（helix）结构，根据茎信息生成图形表示，使用了stem");
draw_helices(stems);

% 初始化缩放功能，允许用户在图形中进行放大或缩小操作
disp("6. 初始化缩放功能，允许用户在图形中进行放大或缩小操作");
setup_zoom();

% 设置图形坐标轴的属性，使整体图形更美观
disp("6. 设置图形坐标轴的属性，使整体图形更美观");
set_nice_axes();

% 更新图形大小，以适应当前窗口或屏幕设置
disp("6. 更新图形大小，以适应当前窗口或屏幕设置");
update_graphics_size();
