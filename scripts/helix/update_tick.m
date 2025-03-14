function update_tick( residue, plot_settings, R )
% update_tick( residue, plot_settings, R )

bp_spacing = plot_settings.bp_spacing;
fontsize = plot_settings.fontsize;

theta = residue.tickrot;
v = [cos(theta*pi/180), sin(theta*pi/180)]*R;
nudge  =  bp_spacing/3;
if length( residue.name ) >= 3 && mod(theta,180) == 90
    nudge = nudge + ( length( residue.name ) - 2 ) * bp_spacing/10;
end
nudge2 =  nudge + bp_spacing/3;
tickpos1 = residue.plot_pos + v*nudge;
tickpos2 = residue.plot_pos + v*nudge2;
set( residue.tick_handle, 'xdata', [tickpos1(1) tickpos2(1)] );
set( residue.tick_handle, 'ydata', [tickpos1(2) tickpos2(2)] );
labelpos = residue.plot_pos + v*nudge2;
set( residue.tick_label, 'position', labelpos );
plot_settings = getappdata( gca, 'plot_settings' );
if ( get(residue.tick_label, 'fontsize') ~= plot_settings.fontsize )
    set( residue.tick_label, 'fontsize', plot_settings.fontsize );
end;
set_text_alignment( residue.tick_label, v );
