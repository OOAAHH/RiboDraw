function plot_pos = get_plot_pos( res_tag, relpos );
if ischar( res_tag ) 
    residue = getappdata( gca, res_tag );
else
    residue = res_tag;
end


%disp('residue fields:');
%disp(fieldnames(residue));
%disp(['residue.chain: ', residue.chain])
%disp(['residue.name: ', residue.name])
%if isfield(residue, 'helix_tag')
    %disp(['residue.helix_tag: ', residue.helix_tag]);
%else
   % error('当前 residue 对象缺少字段 helix_tag。');
%end

helix = getappdata( gca, residue.helix_tag );
R = get_helix_rotation_matrix( helix );
plot_pos = repmat(helix.center,size(relpos,1),1) + relpos*R;
