# RiboDraw 使用教程与命令参考（中文）

这份文档面向「第一次认真用 RiboDraw」的用户：按常见工作流把 **最常用、最容易踩坑** 的命令串起来，并给出可直接复制粘贴的示例。  
RiboDraw 的全部 MATLAB 函数（含内部函数）已经有自动生成的 HTML 文档，但默认不会在 GitHub 里直接展示：请在本地打开 `scripts/docs/menu.html` 浏览完整索引。

---

## 0. 先理解 RiboDraw 的“运行模型”（非常重要）

RiboDraw 绝大多数命令都在 **当前坐标轴 `gca`** 的 `appdata` 里读写数据：

- Residue（残基对象）以 tag 形式保存：`Residue_<chain><segid><resnum>`，例如 `Residue_AA367`
- Helix（螺旋/茎对象）：`Helix_*`
- Linker（连线对象）：`Linker_*`
- Domain/Selection（域/选择）：`Selection_*`
- 全局绘图设置：`plot_settings = getappdata(gca,'plot_settings')`

因此：

- 你必须先有一张 RiboDraw 的 figure/axes（通常由 `initialize_drawing`/`load_drawing` 创建）。
- 若你开了多个 figure，要先切换到正确的图：`figure(n);`（确保命令作用在目标 `gca`）。

常用的“查看当前图里有什么对象”的命令：

```matlab
getappdata(gca)          % 列出当前图 axes 的所有 appdata
get_tags('Residue_')     % 列出所有 Residue 标签
get_tags('Helix_')       % 列出所有 Helix 标签
get_tags('Linker_')      % 列出所有 Linker 标签
gd('Residue_AA367')      % = getappdata(gca,'Residue_AA367')
```

相关函数：`scripts/util/get_tags.m`，`scripts/util/gd.m`。

---

## 1. 安装与路径（Path）设置

### 1.1 把 `scripts/` 加入 MATLAB Path

最推荐的做法是在 RiboDraw 仓库根目录下执行：

```matlab
addpath(genpath(fullfile(pwd,'scripts')));
rehash;
```

检查 MATLAB 找到的函数是否来自你的 RiboDraw 目录：

```matlab
which initialize_drawing -all
which load_drawing -all
which draw_helix -all
```

如果 `which` 显示同名函数来自别的目录（或老版本），说明 path 冲突，需要 `rmpath`/调整顺序。

### 1.2 解释你问到的：`clear ...; rehash;` 到底在干嘛？

当你修改了 `.m` 文件后，MATLAB 有时仍在用“旧版本”的函数（尤其是你把函数文件改坏/修好后）。这时常用：

```matlab
clear draw_helix redraw_helices add_base_pair;
rehash;
```

- `clear funcName`：把该函数从内存里清掉，下次调用会从磁盘重新加载。
- `rehash`：刷新 MATLAB 的文件/路径缓存（新增文件、移动文件、或 `which` 仍指向旧位置时很有用）。

调试时也强烈建议：

```matlab
dbstop if error   % 一旦报错立刻停在出错行，方便定位
```

---

## 2. 输入数据准备：两种启动方式

### 2.1 推荐：用 Rosetta `rna_motif` 生成输入文件

RiboDraw 的标准入口是：

```matlab
initialize_drawing('path/to/your.pdb');
```

其中 `tag = 'path/to/your.pdb'` 只是一个前缀，RiboDraw 会去找同目录下这些文件：

```
your.pdb.fasta
your.pdb.base_pairs.txt
your.pdb.stacks.txt
your.pdb.other_contacts.txt
your.pdb.stems.txt
your.pdb.ligands.txt
your.pdb.motifs.txt
```

它们通常由 Rosetta 的 `rna_motif` 输出得到（详细流程见原始英文教程 `tutorial/tutorial.md`）。

你可以在 `.stems.txt` 的每一行末尾加上 helix 名称（第三列），比如：

```
A:107-110 A:211-214 P4
```

这样螺旋标签会更清晰。

### 2.2 最小化启动：只有序列 + 二级结构（dot-bracket）

如果你没有 3D 结构（也不关心非典型配对/配体/tertiary contacts），至少需要：

- `sequence`：序列字符串
- `secstruct`：dot-bracket 结构字符串
- `resnum_string`：残基编号/链信息，例如 `"A:1-120"`，也支持多段：`"A:1-50 B:1-30"`

先用 `initialize_files` 生成 `.fasta` 和 `.stems.txt`：

```matlab
sequence = 'GGGAAAUCC...';
secstruct = '(((...)))...';
resnum_string = 'A:1-12';
tag = 'myRNA';

initialize_files(sequence, secstruct, resnum_string, tag);
```

然后再运行：

```matlab
initialize_drawing(tag);
```

相关函数：`scripts/setup/initialize_files.m`，`scripts/drawing/initialize_drawing.m`。

---

## 3. 启动与保存：最常用的“生命周期”命令

### 3.1 初始化一张图

```matlab
initialize_drawing('rna_motif/1gid_RNAA.pdb');
```

初始化后常用的“整理视野/交互设置”：

```matlab
set_artboards;         % 视野自动框住所有残基（类似 zoom-to-fit）
show_helix_controls;   % 显示 helix 可拖拽/旋转/反射的控制框
show_linker_controls;  % 显示 linker 顶点控制点（可拖拽调整折线）
```

### 3.2 保存/加载绘图

推荐保存为 `.mat`（快且稳定）：

```matlab
save_drawing('my_drawing.mat');
load_drawing('my_drawing.mat');
```

也支持 JSON（更通用但读写慢）：

```matlab
save_drawing('my_drawing.json');
load_drawing('my_drawing.json');
```

如果你想“把一个 drawing 导入到当前图，但保留当前图里其它没覆盖的对象”，用：

```matlab
load_drawing('subdrawing.mat', 1);   % keep_previous_drawing = 1
```

相关函数：`scripts/drawing/save_drawing.m`，`scripts/drawing/load_drawing.m`。

### 3.3 导出图片（PDF/PNG/SVG 等）

```matlab
export_drawing('my_drawing.pdf');
export_drawing('my_drawing.png');
export_drawing('my_drawing.svg');    % 依赖 plot2svg（若已安装）
```

如果你还想输出每个 residue 的坐标（用于下游软件/对齐），可：

```matlab
export_drawing('my_drawing.png', 1);   % 生成 my_drawing.png.coords.txt
```

相关函数：`scripts/drawing/export_drawing.m`，`scripts/drawing/export_coordinates.m`。

---

## 4. 交互编辑：helix / residue / linker 的常用开关

### 4.1 Helix 控制框（拖动、旋转、反射）

```matlab
show_helix_controls;   % 或 show_helix_controls(1)
hide_helix_controls;   % 等价于 show_helix_controls(0)
```

相关函数：`scripts/helix/show_helix_controls.m`，`scripts/helix/hide_helix_controls.m`，`scripts/helix/draw_helix.m`。

### 4.2 Linker 控制点（折线顶点）

```matlab
show_linker_controls;  % 或 show_linker_controls(1)
hide_linker_controls;  % 等价于 show_linker_controls(0)
```

开启后每条 linker 会出现可拖拽的顶点（小圆点），用于把线段调整为更“横平竖直”的走线。

如果你已经把某两残基之间的 linker 拖出了很多控制点，现在想**还原成默认直线**（清空中间控制点），可以这样做：

```matlab
% 例：还原 A:A:728 与 A:A:729 之间的 arrow linker（按需改成你的残基标签）
res_tag1 = get_res('A:A:770'); res_tag1 = res_tag1{1};
res_tag2 = get_res('A:A:771'); res_tag2 = res_tag2{1};

r1 = getappdata(gca, res_tag1);
r2 = getappdata(gca, res_tag2);
common = intersect(r1.linkers, r2.linkers);

arrow_linkers = get_tags('Linker_','arrow', common);
if isempty(arrow_linkers); error('在这两个残基之间没找到 arrow linker'); end
linker_tag = arrow_linkers{1};

linker = delete_linker(linker_tag, 0); % 只清图形句柄/控制点，不从 drawing 中移除
res1 = getappdata(gca, linker.residue1);
res2 = getappdata(gca, linker.residue2);
linker.relpos1 = res1.relpos;
linker.relpos2 = res2.relpos;
if isfield(linker,'plot_pos'); linker = rmfield(linker,'plot_pos'); end
setappdata(gca, linker.linker_tag, linker);
draw_linker(linker.linker_tag);
```

如果你想让它**自动重新走线**（L 形/折线 heuristic），可以直接：

```matlab
autotrace_linker(linker_tag);
```

如果你想**一键重置整张图的所有 arrow linkers**（清掉所有手工添加的中间控制点，恢复成默认直线），可以：

```matlab
arrow_tags = get_tags('Linker_','arrow');

for i = 1:numel(arrow_tags)
    linker = delete_linker(arrow_tags{i}, 0); % 只清图形句柄/控制点
    res1 = getappdata(gca, linker.residue1);
    res2 = getappdata(gca, linker.residue2);
    linker.relpos1 = res1.relpos;
    linker.relpos2 = res2.relpos;
    if isfield(linker,'plot_pos'); linker = rmfield(linker,'plot_pos'); end
    setappdata(gca, linker.linker_tag, linker);
end

draw_linker(arrow_tags);
```

相关函数：`scripts/settings/show_linker_controls.m`，`scripts/settings/hide_linker_controls.m`。

### 4.3 常见显示开关（show/hide）

这些函数会改 `plot_settings` 并触发重绘（名称非常直观）：

```matlab
show_noncanonical_pairs(1);   % 显示非典型碱基对
show_noncanonical_pairs(0);   % 隐藏

show_stem_pairs(1);           % 显示茎内/典型配对
show_stacks(1);               % 显示 stacking（默认通常关闭）
show_other_contacts(1);       % 显示 other_contacts
show_domains(1);              % 显示 domain 标签
show_motifs(1);               % 显示 motif 相关元素
show_tertiary_contacts(1);    % 显示 tertiary contacts

set_artboards;                % 回到“全图可见”的视野
```

所有开关函数都在：`scripts/settings/`（完整列表见 `scripts/docs/settings/menu.html`）。

---

## 5. 颜色与标注（最常用）

### 5.1 给整张图/某个域/某段残基上色：`color_drawing`

```matlab
color_drawing('rainbow');              % 全图按 resnum 彩虹色（类似 PyMOL spectrum）
color_drawing('black');                % 全图黑色
color_drawing('teal','A:1-50');        % 只给 A:1-50 上色（支持 PyMOL 色名）
color_drawing([1 0 0],'H54b');         % 给某个 helix/selection 上 RGB
```

第二个参数 `selection` 支持：

- `'all'`（默认）
- `'A:1-50'`、`'B:QA:2-10'` 等范围（segid 可选）
- helix 名称（例如 `P4`）或 selection/domain 名称（如果已定义）

相关函数：`scripts/color/color_drawing.m`，`scripts/util/get_res.m`，`scripts/util/get_resnum_from_tag.m`。

---

## 6. 你关心的“手工加一条 H/S/T 的 linker（非典型碱基对）”

你想加的 `A:A:367` 与 `A:A:393` 之间的 `H S T`，对应的是 Leontis-Westhof 表达：

- `edge1 = 'H'`
- `edge2 = 'S'`
- `LW_orientation = 'T'`

最推荐用交互封装好的命令（会自动安装 linker 并触发重绘）：

```matlab
add_base_pair('A:A:367','A:A:393','H','S','T');
```

删除这两者之间的 base-pair 类型 linker：

```matlab
add_base_pair('A:A:367','A:A:393',[],[],[],'delete');
```

相关函数：`scripts/setup/add_base_pair.m`，`scripts/setup/setup_base_pair_linkers.m`。

> 注意  
> `add_base_pair` 依赖 `gca` 里已经有这两个 Residue 对象（也就是你已经 `initialize_drawing` 或 `load_drawing` 过）。

---

## 7. 切片（slice）与合并（merge）：大 RNA 必备

当图很大、拖动很卡时，把某个子域单独开窗口编辑会快很多：

```matlab
slice_drawing('A:120-240');    % 切出一段残基
slice_drawing('Domain I');     % 切出一个 domain（如果你定义过）
```

在子窗口里编辑完后，合并回父图：

```matlab
merge_drawing();               % 默认合并回 slice 时记录的 parent figure
% 或 merge_drawing(1) 指定合并到 figure(1)
```

相关函数：`scripts/drawing/slice_drawing.m`，`scripts/drawing/merge_drawing.m`。

---

## 8. Domain（域）与 Selection（选择）的定义/管理

### 8.1 定义一个 domain（手动指定 residues）

```matlab
setup_domain('A:1-50', 'Domain I');
setup_domain({'A:1-50','A:120-150'}, 'Two blocks');
list_domains;
```

常见操作：

```matlab
show_domains(1);          % 显示 domain 标签
hide_domains;             % 等价 show_domains(0)
show_domain_controls;     % 显示 domain 的可拖拽/编辑控件
hide_domain_controls;
```

相关函数：`scripts/selections/domain/setup_domain.m`，`scripts/selections/domain/list_domains.m`，`scripts/settings/show_domains.m`。

---

## 9. Eterna：导出 customLayout（如果你需要）

把当前图（或某个 selection/domain）导出为可粘贴进 Eterna 的 JSON 片段：

```matlab
export_eterna();              % 全图
export_eterna('Domain I');    % 只导出某个域
```

更复杂的用法（全长序列/locks/IUPAC/结构约束）见函数头注释：
`scripts/eterna/export_eterna.m`，以及英文教程 `tutorial/eterna/eterna_tutorial.md`。

---

## 10. 常见问题排查（建议收藏）

### 10.1 “Undefined function …” / 找不到命令

1) 先检查路径：

```matlab
which initialize_drawing -all
```

2) 确保加了 `scripts/` 以及所有子目录：

```matlab
addpath(genpath(fullfile(pwd,'scripts')));
rehash;
```

### 10.2 你刚修改了某个 `.m`，但 MATLAB 还在跑旧版本

```matlab
clear functionName;
rehash;
```

典型例子（你刚遇到的）：

```matlab
clear draw_helix redraw_helices add_base_pair;
rehash;
```

### 10.3 报错“函数定义放置或嵌套错误 / 运算符使用无效”

这类错误通常意味着某个 `.m` 文件被合并/编辑坏了（括号/`end` 不匹配、插入了乱码、或有残留冲突标记）。

排查步骤：

```matlab
which draw_helix -all      % 确认你调用的是哪一个文件
edit draw_helix            % 打开检查是否有 <<<<<<< / ======= / >>>>>>>
```

你也可以在仓库里搜索冲突标记：

```bash
rg -n "<<<<<<<|=======|>>>>>>>" scripts
```

---

## 11. 全量函数索引在哪里？

- **按目录分类的 HTML 文档**：`scripts/docs/menu.html`
  - 每个子模块也有自己的菜单页，例如：
    - `scripts/docs/drawing/menu.html`
    - `scripts/docs/settings/menu.html`
    - `scripts/docs/linkers/menu.html`
    - `scripts/docs/helix/menu.html`
- **在 MATLAB 里看函数用法**：
  - `help functionName`
  - `type functionName`（直接看源码）

---

## 12. 命令速查表（按模块）

这一节是“查名字用”的速查表：每个命令后面都给一个最常见的调用方式。更完整的参数说明请看对应 `.m` 文件头注释或 `scripts/docs/` HTML。

### 12.1 Drawing / 读写与导出（`scripts/drawing/`）

- 初始化：`initialize_drawing(tag)`
- 保存：`save_drawing('drawing.mat')` / `save_drawing('drawing.json')`
- 加载：`load_drawing('drawing.mat')`
- 导入（保留当前图的其它对象）：`import_drawing('sub.mat')`
- 导出图片：`export_drawing('out.pdf')` / `export_drawing('out.png', 1)`
- 导出坐标（一般用 `export_drawing(...,1)` 即可）：`export_coordinates('out.coords.txt')`
- 切片：`slice_drawing('A:1-100')` / `slice_drawing('Domain I')`
- 合并：`merge_drawing()` / `merge_drawing(figureNumber)`
- 取二级结构：`get_secstruct_from_drawing()` / `get_secstruct_from_drawing('Domain I')`
- 取序列：`get_sequence_from_drawing()`

### 12.2 Setup / 从文本信息建对象（`scripts/setup/`）

- 只用 sequence/secstruct 启动前的文件生成：`initialize_files(sequence, secstruct, resnum_string, tag)`
- 安装碱基对 linkers：`setup_base_pair_linkers(base_pairs)`
- 安装箭头 linkers：`setup_arrow_linkers(resnum, chains, segid)`
- 安装 stacking/other_contacts：`setup_base_stack_linkers(base_stacks)`、`setup_other_contact_linkers(other_contacts)`
- 交互式补一对碱基对：`add_base_pair(tag1, tag2, edge1, edge2, LW, mode)`

### 12.3 Helix / 螺旋绘制与重绘（`scripts/helix/`）

- 画全部 helices：`draw_helices(stems)`（通常初始化/加载时已做）
- 重绘全部 helices：`redraw_helices()`
- 重绘某个 helix：`redraw_helix(h)`（一般由交互控件触发）
- 显示/隐藏 helix 控制框：`show_helix_controls(1/0)`、`hide_helix_controls`
- 显示/隐藏 helix 标签：`show_helix_label`、`hide_helix_label`

### 12.4 Linkers / 连线（`scripts/linkers/` + `scripts/settings/`）

- 手工把 linker 关联到 residue（通常不需要直接用）：`add_linker_to_residue(res_tag, linker_tag)`
- 删除 linker：`delete_linker(linker)`（可传 tag 或 struct 或 cell）
- 重绘/重算一条 linker：`draw_linker(linker_tag)`、`rectify_linker(linker_tag)`
- 自动走线：`autotrace_linker(linker_tag)`（以及其它 `autotrace_*`）
- 显示/隐藏 linker 顶点控制：`show_linker_controls(1/0)`、`hide_linker_controls`
- 显示/隐藏类别：`show_stem_pairs`、`show_noncanonical_pairs`、`show_stacks`、`show_other_contacts`、`show_tertiary_contacts`

### 12.5 Settings / 显示开关与主题（`scripts/settings/`）

- 视野：`set_artboards`
- 主题：`set_default_theme`
- 背景/线色/符号色：`set_bg_color('white')`、`set_line_color('k')`、`set_symbol_color('k')`
- 字体/线宽：`set_fontsize(12)`、`set_linker_width(fontsize)`
- ticks：`set_tick_frequency(10)`、`set_chain_ticks(1)`、`show_chain_termini(1/0)`
- 图像显示：`show_images(1/0)`、`show_images_as_boundaries`、`show_images_as_rounded_rectangles`

### 12.6 Selections & Domains（`scripts/selections/`）

- 定义 domain：`setup_domain('A:1-50', 'Domain I')`
- 列出 domain：`list_domains`
- 显示/隐藏 domain：`show_domains(1/0)`、`hide_domains`
- 编辑控件：`show_domain_controls`、`hide_domain_controls`、`show_selection_controls`、`hide_selection_controls`

### 12.7 Motifs / Ligands / Tertiary contacts

- Motifs：`setup_motifs(motifs)`、`draw_motifs(motif_tags)`
- Ligands：`show_ligand_linkers`、`hide_ligand_linkers`、`set_ligand_image_scaling`
- Tertiary contacts：`setup_tertiary_contacts`、`clear_tertiary_contacts`、`delete_tertiary_contact`

### 12.8 GUI（`scripts/gui/`）

- `setup_zoom()`、`setup_pan()`（一般初始化/加载时已设置）
- `set_nice_axes()`、`update_graphics_size()`、`move_stuff_to_back()`

### 12.9 Util（`scripts/util/`）

- 解析残基范围字符串：`get_resnum_from_tag('A:1-10')`
- 把选择（resrange/domain/helix 名）转成 res_tags：`get_res('A:1-10')`
- 列出对象 tags：`get_tags('Residue_')` / `get_tags('Linker_')`
- 读对象：`gd(tag)`（= `getappdata(gca,tag)`）

---

如果你愿意，我也可以按你实验室的常见任务（比如“加/删非典型配对、定义 domain、导出 Eterna customLayout、处理 pseudoknot、批量导出图片”）把这份教程再拆成几个更短的“操作手册”，每个手册配一段可复用的脚本模板。
