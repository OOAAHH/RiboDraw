function draw_helices( helices )
% draw_helices( helices )
%
% 渲染 RiboDraw 绘图中初始的螺旋（helix）结构，同时绘制相关的残基、连接器和碱基对。
%
% 该函数是 draw_helix() 的封装（wrapper）。
%
% 如果未提供 helices 参数，则在当前坐标轴（gca）中查找名称包含 "Helix_" 的对象，并将其作为螺旋对象进行绘制。
%
% (C) R. Das, Stanford University, 2017

%% 如果未传入 helices 参数，则从 gca 的应用数据中提取所有螺旋对象
if ~exist( 'helices', 'var' )
    fprintf('helices 参数未传入，开始从 gca 中提取螺旋对象...\n');
    helices = {};                         % 初始化空 cell 数组
    appdata = getappdata( gca );            % 获取当前坐标轴上所有的应用数据
    all_fields = fields(appdata);
    fprintf('gca 中共有 %d 个应用数据字段\n', length(all_fields));
    objnames = all_fields;                % 获取应用数据中的所有字段名称
    for n = 1:length( objnames )
        % 如果字段名称中包含 "Helix_" 字符串，则认为该字段对应一个螺旋对象
        if ~isempty( strfind( objnames{n}, 'Helix_' ) );
            % 使用 getfield 从 appdata 中取出该螺旋对象，并加入 helices 数组
            helices = [helices, getfield(appdata, objnames{n})];
            fprintf('提取到螺旋对象：%s\n', objnames{n});
        end
    end
    fprintf('共提取到 %d 个螺旋对象。\n', length(helices));
else
    fprintf('helices 参数已传入，共 %d 个对象。\n', length(helices));
end

%% 临时隐藏域（domains）的显示以提高绘图速度
% 保存当前的显示状态，并将其临时关闭，待绘制完后再恢复显示
fprintf('开始临时隐藏域显示...\n');
save_show_domains = temporarily_hide_domains();
fprintf('临时隐藏前，域显示状态为：%d\n', save_show_domains);

%% 绘制每个螺旋结构
% 使用 textprogressbar 显示绘制进度信息
fprintf('开始绘制螺旋结构...\n');
textprogressbar('Drawing helices... ');
for n = 1:length( helices )
    fprintf('正在绘制第 %d 个螺旋...\n', n);
    % 对每个螺旋调用 draw_helix 进行绘制
    draw_helix( helices{n} );
    % 更新进度条显示当前绘制进度（百分比）
    textprogressbar( 100 * n/length(helices) );
end
textprogressbar(' done');  % 绘制结束，更新进度条状态
fprintf('所有螺旋结构绘制完成。\n');
axis off;                   % 关闭坐标轴显示

%% 设置图形背景颜色
% 默认背景颜色为白色
bg_color = 'white';
plot_settings = get_plot_settings();
% 如果绘图设置中指定了背景颜色，则使用该颜色
if isfield( plot_settings, 'bg_color' )
    bg_color = plot_settings.bg_color;
end;
fprintf('设置图形背景颜色为：%s\n', bg_color);
set(gcf, 'color', bg_color);  % 设置当前图形窗口的背景颜色

%% 恢复域（domains）的显示
fprintf('检查是否需要恢复域显示...\n');
if save_show_domains
    fprintf('恢复域显示...\n');
    show_domains(1, 0);  % 恢复域的显示
end;

%% 如果设置中要求显示“base rope”，则调用相应函数绘制
if isfield( plot_settings, 'show_base_rope' ) && plot_settings.show_base_rope
    fprintf('绘制 base rope...\n');
    draw_base_rope();
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 嵌套函数：暂时隐藏域的显示
function save_show_domains = temporarily_hide_domains();
    % 从当前坐标轴的应用数据中获取绘图设置
    plot_settings = getappdata( gca, 'plot_settings' );
    % 保存当前域显示的状态
    save_show_domains = plot_settings.show_domains;
    fprintf('当前域显示状态：%d\n', save_show_domains);
    % 如果当前设置中域是显示的，则将其关闭，并更新绘图设置
    if save_show_domains
        fprintf('关闭域显示以提高绘图速度...\n');
        plot_settings.show_domains = 0;
        setappdata( gca, 'plot_settings', plot_settings );
    end;
