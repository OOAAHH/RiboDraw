function show_chain_termini( setting )
% show_chain_termini( setting )
%
% Show or hide 5'/3' labels at chain termini.
%
% (C) R. Das, Stanford University, 2017

if ~exist( 'setting', 'var' ); setting = 1; end;

plot_settings = getappdata( gca, 'plot_settings' );
plot_settings.show_chain_termini = setting;
setappdata( gca, 'plot_settings', plot_settings );

draw_helices();

