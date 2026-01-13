function helix = draw_helix( helix )
% helix = draw_helix( helix )
%
% RiboDraw 的主绘图函数，负责绘制一个“螺旋”结构，
% 包括由 Watson-Crick 配对残基组成的茎，以及附近的环状残基，
% 这些残基会随着螺旋一起平移、反射或旋转。
%
% 此函数还更新了与这些残基相关联的连接器、选择（域）、螺旋标签和刻度。
%
% 注意：名称 “helix” 可能并不完全准确，源于最初的草稿代码。
%
% TODO: 更新以处理非螺旋型 RNA 结构（例如 G-quadruplexes）。
%
% 输入：
%   helix = 螺旋对象，或当前绘图（gca）中螺旋对象的标签（字符串）
%
% 输出：
%   helix = 更新后的螺旋对象
%
% (C) R. Das, Stanford University, 2017-2018

% 如果输入的 helix 为字符串，则从 gca 中获取相应的对象
if ischar( helix )
    fprintf('【调试】输入 helix 为字符串，尝试通过标签获取对象：%s\n', helix);
    helix = getappdata(gca, helix);
end

% 获取全局绘图设置
plot_settings = getappdata( gca, 'plot_settings' );
fprintf('【调试】绘图设置已获取。\n');

% 提取螺旋中心和旋转矩阵
helix_center = helix.center;
R = get_helix_rotation_matrix( helix ); 
fprintf('【调试】螺旋中心： [%s]\n', num2str(helix_center));
fprintf('【调试】旋转矩阵 R：\n');
disp(R);

% 获取螺旋茎中第一个配对残基的数量
N = length( helix.resnum1 );
spacing = plot_settings.spacing;
bp_spacing = plot_settings.bp_spacing;
helix_res_tags = {};

% 循环处理螺旋茎上每对配对残基，计算它们的相对位置
for k = 1:N
    % 第一个配对残基（下侧），构造标签并更新位置
    res_tag = sanitize_tag(sprintf('Residue_%s%s%d', helix.chain1(k), helix.segid1{k}, helix.resnum1(k)));
    fprintf('【调试】处理配对残基1，第 %d 个，标签：%s\n', k, res_tag);
    pos1 = update_residue_pos( res_tag, [ spacing*((k-1)-(N-1)/2), -bp_spacing/2], helix.center, R );
    helix_res_tags = [helix_res_tags, res_tag ];
    
    % 第二个配对残基（上侧），注意逆序处理
    res_tag = sanitize_tag(sprintf('Residue_%s%s%d', helix.chain2(N-k+1), helix.segid1{N-k+1}, helix.resnum2(N-k+1)));
    fprintf('【调试】处理配对残基2，第 %d 个，标签：%s\n', k, res_tag);
    pos2 = update_residue_pos( res_tag, [ spacing*((k-1)-(N-1)/2), +bp_spacing/2], helix.center, R );
    helix_res_tags = [helix_res_tags, res_tag ];
    
    all_pos1(k,:) = pos1;
    all_pos2(k,:) = pos2;
    
    fprintf('【调试】第 %d 对配对残基位置：pos1 = [%s]，pos2 = [%s]\n', k, num2str(pos1), num2str(pos2));
end

% 绘制与该螺旋相关联的所有残基
not_helix_res_tags = {};
fprintf('【调试】开始绘制与螺旋关联的残基...\n');
for i = 1:length( helix.associated_residues )
    res_tag = helix.associated_residues{i};
    fprintf('【调试】处理关联残基：%s\n', res_tag);
    residue = getappdata( gca, res_tag );
    if ~isfield( residue, 'name' )
        fprintf('【警告】残基 %s 不包含字段 name，跳过。\n', res_tag);
        continue;
    end
    if ~isfield( residue, 'relpos' )
        fprintf('【调试】残基 %s 缺少 relpos，设置默认相对位置。\n', res_tag);
        residue.relpos = set_default_relpos( residue, helix, plot_settings );
        setappdata( gca, res_tag, residue );
    end
    draw_residue_for_helix( res_tag, helix_center, R, plot_settings );
    % 如果该残基不在螺旋茎内，则加入 not_helix_res_tags
    if ~any(strcmp(helix_res_tags, res_tag))
        not_helix_res_tags = [not_helix_res_tags, res_tag];
    end;
end

% 更新与这些残基相关联的所有连接器（linkers）
redrawn_linkers = {};
fprintf('【调试】开始更新连接器...\n');

% 在开始更新时，打印当前 gca 中所有 linker 的名称
appdata = getappdata(gca);
all_fields = fieldnames(appdata);
fprintf('【调试】当前所有 linker 名称：\n');
for j = 1:length(all_fields)
    if ~isempty(strfind(all_fields{j}, 'Linker_'))
        fprintf('  %s\n', all_fields{j});
    end
end

for i = 1:length(helix.associated_residues)
    res_tag = helix.associated_residues{i};
    residue = getappdata(gca, res_tag);
    
    % 打印当前残基所有字段信息
    fprintf('【调试】当前残基 %s 内部所有字段：\n', res_tag);
    disp(residue);
    
    if ~isfield(residue, 'linkers')
        fprintf('【调试】残基 %s 无连接器字段，跳过。\n', res_tag);
        continue;
    end
    linker_tags = residue.linkers;
    for k = 1:length(linker_tags)
        if any(strcmp(redrawn_linkers, linker_tags{k}))
            fprintf('【调试】连接器 %s 已更新，跳过重复绘制。\n', linker_tags{k});
            continue;
        end
        linker = getappdata(gca, linker_tags{k});
        
        % 打印当前链接器内部所有字段信息
        fprintf('【调试】当前链接器 %s 内部所有字段：\n', linker_tags{k});
        disp(linker);
        
        % 再次打印当前残基信息，便于确认更新前状态
        fprintf('【调试】更新前，残基 %s 内部所有字段：\n', res_tag);
        disp(residue);
        
        fprintf('【调试】更新连接器：%s\n', linker_tags{k});
        draw_linker(linker, plot_settings);
        redrawn_linkers = [redrawn_linkers, linker.linker_tag];
    end
end




% 如果 helix 中未包含 label_relpos 字段，则使用 bp_spacing 设置默认值
if ~isfield( helix, 'label_relpos' )
    helix.label_relpos = plot_settings.bp_spacing * [0 1];
    fprintf('【调试】设置默认 label_relpos：[%s]\n', num2str(helix.label_relpos));
end
% 更新螺旋标签，并返回更新后的 helix 对象
helix = make_helix_label( helix, plot_settings, R );
fprintf('【调试】螺旋标签已更新。\n');

% 处理 selections 和 motifs（域与模体）的绘制
selections = {};
motifs = {};
for i = 1:length( helix.associated_residues )
    res_tag = helix.associated_residues{i};
    residue = getappdata( gca, res_tag );
    if isfield( residue, 'associated_selections' ) && ~isempty( residue.associated_selections )
        selections = [ selections, residue.associated_selections ];
    end    
    if isfield( residue, 'associated_motifs' ) && ~isempty( residue.associated_motifs )
        motifs = [ motifs, residue.associated_motifs ];
    end    
end
selections = unique( selections );
fprintf('【调试】共检测到 %d 个 selections。\n', length(selections));
draw_selections( selections );
motifs = unique( motifs );
fprintf('【调试】共检测到 %d 个 motifs。\n', length(motifs));
draw_motifs( motifs );

% 为螺旋编辑创建图形句柄 —— 矩形（用于拖动）
minpos = min( [all_pos1; all_pos2] );
maxpos = max( [all_pos1; all_pos2] );
fprintf('【调试】计算螺旋矩形范围：minpos = [%s]，maxpos = [%s]\n', num2str(minpos), num2str(maxpos));
helix = create_default_rectangle( helix, 'helix_tag', helix.helix_tag, @redraw_helix );
set_rectangle_coords( helix, minpos, maxpos, spacing );

% 为螺旋添加反射用的可点击线条
if ~isfield( helix, 'reflect_line1' )
    h = plot( [0 0], [0 0], 'color',[0.5 0.5 1],'clipping','off' );
    setappdata( h, 'helix_tag', helix.helix_tag);
    set(h,'ButtonDownFcn',{@reflect_helix,h});
    helix.reflect_line1 = h;
    fprintf('【调试】创建反射线1。\n');
end
line1 = helix_center + spacing*[-(N+0.25)/2, 0]*R;
line1x = helix_center + spacing*[-(N-0.75)/2, 0]*R;
set( helix.reflect_line1, 'Xdata', [line1(1) line1x(1)], 'Ydata', [line1(2) line1x(2)]);
fprintf('【调试】设置反射线1位置。\n');

if ~isfield( helix, 'reflect_line2' )
    h = plot( [0 0], [0 0], 'color',[0.5 0.5 1],'clipping','off' );
    setappdata( h, 'helix_tag', helix.helix_tag);
    set(h,'ButtonDownFcn',{@reflect_helix,h});
    helix.reflect_line2 = h;
    fprintf('【调试】创建反射线2。\n');
end
line2 = helix_center + spacing*[ (N+0.25)/2, 0]*R;
line2x = helix_center + spacing*[ (N-0.75)/2, 0]*R;
set( helix.reflect_line2, 'Xdata', [line2(1) line2x(1)], 'Ydata', [line2(2) line2x(2)]);
fprintf('【调试】设置反射线2位置。\n');

% 为螺旋添加可点击的旋转中心
if ~isfield( helix, 'click_center' )
    h = rectangle( 'Position',[0 0 0 0], ...
                   'curvature',[0.5 0.5],...
                   'edgecolor',[0.5 0.5 1],...
                   'facecolor',[0.5 0.5 1],...
                   'linewidth',1.5,'clipping','off' );
    setappdata( h,'helix_tag', helix.helix_tag);
    set(h,'ButtonDownFcn',{@rotate_helix,h});
    helix.click_center = h;
    fprintf('【调试】创建旋转中心对象。\n');
end
set( helix.click_center, 'Position', [helix_center(1)-0.15*spacing, helix_center(2)-0.15*spacing, 0.3*spacing, 0.3*spacing]);
fprintf('【调试】设置旋转中心位置。\n');

% 设置螺旋控制界面是否可见
set_helix_visibility( helix, plot_settings.show_helix_controls );

% 设置刻度标签为可拖动
for i = 1:length( helix.associated_residues )
    res_tag = helix.associated_residues{i};
    residue = getappdata( gca, res_tag );
    if isfield( residue, 'tick_label' ) && isvalid( residue.tick_label )
        setappdata( residue.tick_label, 'res_tag', res_tag );
        if ~isappdata(residue.tick_label,'user_movefcn')
            draggable( residue.tick_label, @move_tick, 'endfcn', @redraw_tick_res_and_helix );
            fprintf('【调试】设置刻度标签 %s 为可拖动。\n', res_tag);
        end;
    end
end

% 设置单链残基为可拖动
for i = 1:length( not_helix_res_tags )
    res_tag = not_helix_res_tags{i};
    residue = getappdata( gca, res_tag );
    if ~isappdata(residue.handle,'user_movefcn')
        draggable( residue.handle, 'n', [-inf inf -inf inf], @move_snapgrid, 'endfcn', @redraw_res_and_helix );
        fprintf('【调试】设置非螺旋残基 %s 为可拖动。\n', res_tag);
    end;
end

% 绘制基础连接线（base rope）
draw_base_rope();
fprintf('【调试】绘制 base rope 完成。\n');

%%%%%%%%%%%%%%%%%%%%%
% 最后将更新后的 helix 对象存回 gca 中
setappdata( gca, helix.helix_tag, helix );
fprintf('【调试】更新后的 helix 对象已存回 gca，标签：%s\n', helix.helix_tag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下为内部辅助函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：绘制与螺旋相关的残基
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = draw_residue_for_helix( res_tag, helix_center, R, plot_settings )
    residue = getappdata( gca, res_tag );
    if isfield( residue, 'relpos' )
        pos = helix_center + residue.relpos * R;
        residue.plot_pos = pos;
        residue.res_tag = res_tag;
        fprintf('【调试】绘制残基 %s, 计算位置：[%s]\n', res_tag, num2str(pos));
        residue = draw_residue( residue );
        h = residue.handle;
        set( h, 'Position', pos );
        setappdata( residue.handle, 'res_tag', res_tag );
        residue = draw_tick( residue, plot_settings, R );
        residue = draw_chain_termini_labels( residue, plot_settings, R );
        setappdata( gca, res_tag, residue );
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：绘制链起止(5'/3')标签
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residue = draw_chain_termini_labels( residue, plot_settings, R )
    if ~isfield( plot_settings, 'show_chain_termini' ) || ~plot_settings.show_chain_termini
        residue = rmgraphics( residue, {'five_prime_label','three_prime_label','five_prime_tick_label','five_prime_tick_handle','three_prime_tick_label','three_prime_tick_handle'} );
        return;
    end

    if ~isfield( residue, 'plot_pos' ) || isempty( residue.plot_pos )
        return;
    end

    color = 'k';
    if isfield( plot_settings, 'line_color' )
        color = plot_settings.line_color;
    end

    % Use tick-style positioning so labels move like residue number ticks.
    if ~isfield( residue, 'tickrot' ) || isempty( residue.tickrot ) || isnan( residue.tickrot )
        residue = set_default_tickrot( residue );
    end

    show_five_prime  = isfield( residue, 'is_five_prime'  ) && residue.is_five_prime;
    show_three_prime = isfield( residue, 'is_three_prime' ) && residue.is_three_prime;

    five_theta = residue.tickrot;
    three_theta = residue.tickrot;
    if show_five_prime && show_three_prime
        three_theta = mod( residue.tickrot + 180, 360 );
    end

    if show_five_prime
        [residue.five_prime_tick_handle, residue.five_prime_tick_label] = draw_chain_terminus_tick( ...
            residue, R, plot_settings, color, five_theta, 'five_prime_tick_handle', 'five_prime_tick_label', '5''' );
    else
        residue = rmgraphics( residue, {'five_prime_tick_label','five_prime_tick_handle'} );
    end

    if show_three_prime
        [residue.three_prime_tick_handle, residue.three_prime_tick_label] = draw_chain_terminus_tick( ...
            residue, R, plot_settings, color, three_theta, 'three_prime_tick_handle', 'three_prime_tick_label', '3''' );
    else
        residue = rmgraphics( residue, {'three_prime_tick_label','three_prime_tick_handle'} );
    end

    % Legacy cleanup (older label field names).
    residue = rmgraphics( residue, {'five_prime_label','three_prime_label'} );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：绘制一个末端 tick（细线 + 文本）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tick_handle, tick_label] = draw_chain_terminus_tick( residue, R, plot_settings, color, theta, handle_field, label_field, label_string )
    if isfield( residue, handle_field ) && isvalid( getfield( residue, handle_field ) )
        tick_handle = getfield( residue, handle_field );
    else
        tick_handle = plot( [0, 0], [0, 0], color, 'linewidth', 0.5, 'clipping', 'off' );
    end

    if isfield( residue, label_field ) && isvalid( getfield( residue, label_field ) )
        tick_label = getfield( residue, label_field );
    else
        tick_label = text( 0, 0, label_string, ...
            'fontsize', plot_settings.fontsize, ...
            'horizontalalign', 'center', 'verticalalign', 'middle', ...
            'clipping', 'off', 'color', color );
    end

    v = [cos(theta*pi/180), sin(theta*pi/180)] * R;
    bp_spacing = plot_settings.bp_spacing;
    nudge = bp_spacing/3;
    if length( residue.name ) >= 3 && mod(theta,180) == 90
        nudge = nudge + ( length( residue.name ) - 2 ) * bp_spacing/10;
    end
    nudge2 = nudge + bp_spacing/3;
    tickpos1 = residue.plot_pos + v*nudge;
    tickpos2 = residue.plot_pos + v*nudge2;
    set( tick_handle, 'xdata', [tickpos1(1) tickpos2(1)] );
    set( tick_handle, 'ydata', [tickpos1(2) tickpos2(2)] );

    labelpos = residue.plot_pos + v*nudge2;
    set( tick_label, 'position', labelpos );
    if ( get(tick_label, 'fontsize') ~= plot_settings.fontsize )
        set( tick_label, 'fontsize', plot_settings.fontsize );
    end
    set( tick_label, 'String', label_string );
    set_text_alignment( tick_label, v );

    setappdata( tick_label, 'res_tag', residue.res_tag );
    if ~isappdata( tick_label, 'user_movefcn' )
        draggable( tick_label, @move_tick, 'endfcn', @redraw_tick_res_and_helix );
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：设置默认相对位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function relpos = set_default_relpos( residue, helix, plot_settings )
    dist1 = min( abs(helix.resnum1 - residue.resnum) );
    if ( residue.chain ~= helix.chain1(1) )
        dist1 = Inf * dist1;
    end
    dist2 = min( abs(helix.resnum2 - residue.resnum) );
    if ( residue.chain ~= helix.chain2(1) )
        dist2 = Inf * dist2;
    end
    [~, strand] = min( [min(dist1), min(dist2)] );
    N = length( helix.resnum1 );
    if ( strand == 1 )
        d = residue.resnum - helix.resnum1(1);
        if abs(d) > 10
            d = sign(d) * 10 * ( log( abs(d)/ 10) + 1 );
        end
        relpos = [ plot_settings.spacing*(d - (N-1)/2), -plot_settings.bp_spacing/2 ];
    else
        assert( strand == 2 );
        d = residue.resnum - helix.resnum2(1);
        if abs(d) > 10
            d = sign(d) * 10 * ( log( abs(d)/ 10) + 1 );
        end
        relpos = [ plot_settings.spacing*(-d + (N-1)/2), +plot_settings.bp_spacing/2 ];
    end
    fprintf('【调试】设置默认 relpos 为：[%s]\n', num2str(relpos));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：更新残基位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = update_residue_pos( res_tag, relpos, center, R )
    residue = getappdata( gca, res_tag );
    residue.relpos = relpos;
    pos = center + relpos * R;
    setappdata( gca, res_tag, residue);
    fprintf('【调试】更新残基 %s 的 relpos 为：[%s]，计算位置 pos = [%s]\n', res_tag, num2str(relpos), num2str(pos));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：生成螺旋标签
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function helix = make_helix_label( helix, plot_settings, R )
    % 如果不存在 label 或 label 对象无效，则创建新的文本对象
    if ~isfield( helix, 'label' ) || ~isvalid( helix.label)
        h = text( 0, 0, helix.name, 'fontsize', plot_settings.fontsize*1.5, ...
                  'fontname','helvetica','clipping','off');
        helix.label = h;
        setappdata( helix.label, 'helix_tag', helix.helix_tag );
        draggable( helix.label, 'n', [-inf inf -inf inf], @move_helix_label, 'endfcn', @redraw_helix_label );
        fprintf('【调试】创建螺旋标签对象。\n');
    end
    h = helix.label;
    label_pos = helix.center + helix.label_relpos * R;
    set( h, 'String', helix.name );
    set( h, 'position', label_pos );
    set( h, 'fontsize', plot_settings.fontsize*1.5 );
    color = 'k';
    if isfield( plot_settings, 'line_color' )
        color = plot_settings.line_color;
    end
    if isfield( helix, 'rgb_color' )
        color = helix.rgb_color;
    end
    set( h, 'color', color );
    v = [0, sign(helix.label_relpos(2))] * R;
    set_text_alignment( h, v );
    if isfield( helix, 'label_visible' )
        if helix.label_visible
            visible = 'on';
        else
            visible = 'off';
        end
        set( helix.label, 'visible', visible );
    end
    fprintf('【调试】螺旋标签位置设置为：[%s]\n', num2str(label_pos));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：移动螺旋标签
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function move_helix_label(h)
    pos = get(h, 'position'); 
    helix_tag = getappdata( h, 'helix_tag' );
    helix = getappdata( gca, helix_tag );
    R = get_helix_rotation_matrix( helix );
    relpos = (pos(1:2) - helix.center) * R';
    v = [0, sign(relpos(2))] * R;
    set_text_alignment( h, v );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：重绘螺旋标签
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_helix_label(h)
    pos = get(h, 'position'); 
    helix_tag = getappdata( h, 'helix_tag' );
    helix = getappdata( gca, helix_tag );
    R = get_helix_rotation_matrix( helix );
    helix.label_relpos = (pos(1:2) - helix.center) * R';
    plot_settings = getappdata( gca, 'plot_settings' );
    snap_spacing = plot_settings.bp_spacing/4;
    helix.label_relpos = round( helix.label_relpos / snap_spacing ) * snap_spacing;
    draw_helix( helix );
    fprintf('【调试】重绘螺旋标签完成。\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：绘制刻度（tick）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residue = draw_tick( residue, plot_settings, R )
    if mod(residue.resnum, plot_settings.tick_frequency) ~= 0
        % 清除旧刻度信息
        if isfield( residue, 'tickrot' )
            residue = rmfield( residue, 'tickrot' );
        end
        if isfield( residue, 'tick_handle' )
            residue = rmfield( residue, 'tick_handle' );
        end
        if isfield( residue, 'tick_label' )
            residue = rmfield( residue, 'tick_label' );
        end
        return;
    end
    if isfield(residue, 'ligand_partners')
        return;
    end
    if ~isfield( residue, 'tickrot' )
        residue.tickrot = nan; % 稍后根据螺旋旋转设置
    end
    color = 'k';
    if isfield( plot_settings, 'line_color' )
        color = plot_settings.line_color;
    end
    if ~isfield( residue, 'tick_handle' ) || ~isvalid( residue.tick_handle )
        residue.tick_handle = plot( [0, 0], [0, 0], color, 'linewidth', 0.5, 'clipping', 'off');
        setappdata( gca, residue.res_tag, residue );
    end
    if plot_settings.chain_ticks == 1
        residue.tick_label = text( 0, 0, sprintf('%d(%s)', residue.resnum, residue.chain), ...
            'fontsize', plot_settings.fontsize, 'horizontalalign', 'center', ...
            'verticalalign', 'middle', 'clipping', 'off', 'color', color );
        setappdata( gca, residue.res_tag, residue );
    end
    if ~isfield( residue, 'tick_label' ) || ~isvalid( residue.tick_label )
        if plot_settings.chain_ticks == 1
            residue.tick_label = text( 0, 0, sprintf('%d(%s)', residue.resnum, residue.chain), ...
                'fontsize', plot_settings.fontsize, 'horizontalalign', 'center', ...
                'verticalalign', 'middle', 'clipping', 'off', 'color', color );
        else
            residue.tick_label = text( 0, 0, num2str(residue.resnum), ...
                'fontsize', plot_settings.fontsize, 'horizontalalign', 'center', ...
                'verticalalign', 'middle', 'clipping', 'off', 'color', color );
        end
        setappdata( gca, residue.res_tag, residue );
    end
    if isfield( residue, 'tickrot' )
        if isnan(residue.tickrot)
            residue = set_default_tickrot( residue );
        end
        update_tick( residue, plot_settings, R );
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：设置默认的刻度旋转角度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residue = set_default_tickrot( residue )
    if sign( residue.relpos(2) ) > 0 
        residue.tickrot = 90;
    else
        residue.tickrot = 270;
    end
    fprintf('【调试】设置残基 %s 的默认 tickrot 为 %d\n', residue.res_tag, residue.tickrot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：拖动刻度时的处理（snap to grid）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function move_tick(h)
    pos = get(h, 'Position');
    pos = pos(1:2);
    residue = getappdata( gca, getappdata(h, 'res_tag') );
    v = pos - residue.plot_pos;
    theta = atan2( v(2), v(1) );
    theta = round(theta * 180/pi/45) * 45;
    v = [cos(theta*pi/180), sin(theta*pi/180)];
    plot_settings = getappdata(gca, 'plot_settings');
    labelpos = residue.plot_pos + v * plot_settings.bp_spacing * 2/3;
    set(h, 'Position', labelpos);
    set_text_alignment( h, v );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：重绘刻度和螺旋
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_tick_res_and_helix(h)
    pos = get(h, 'Position');
    pos = pos(1:2);
    res_tag = getappdata(h, 'res_tag');
    residue = getappdata( gca, res_tag );
    v = pos - residue.plot_pos;
    helix = getappdata( gca, residue.helix_tag );
    R = get_helix_rotation_matrix( helix );
    v = v * R';
    tickrot = mod(round(atan2(v(2), v(1))*180/pi/45)*45, 360);
    residue.tickrot = tickrot;
    setappdata( gca, res_tag, residue );
    redraw_res_and_helix( residue.handle );
    fprintf('【调试】重绘刻度和螺旋完成，res_tag：%s\n', res_tag);

% 以上为 draw_helix 及其辅助函数的全部代码
