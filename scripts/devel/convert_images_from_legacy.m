function convert_images_from_legacy()

plot_settings = getappdata( gca, 'plot_settings' );
res_tags = get_tags( 'Residue' );
for i = 1:length( res_tags )
    ligand = getappdata( gca, res_tags{i} );
    if ~isfield( ligand, 'ligand_partners' ); continue; end
    %ligand = draw_image_boundary( ligand, plot_settings );
    if( isfield( ligand, 'image_offset' ) & ~all(ligand.image_offset == 0 ) )
        ligand.plot_pos      = get_plot_pos( ligand, ligand.relpos ) + ligand.image_offset;
        ligand.relpos        = get_relpos( ligand.plot_pos, getappdata(gca,ligand.helix_tag) );
        ligand.label_relpos  = get_relpos( ligand.plot_pos - ligand.image_offset, getappdata(gca,ligand.helix_tag) ); 
        ligand.image_offset  = [0,0];
        ligand = rmfield( ligand, 'image_offset' );
    end    
    setappdata(gca,res_tags{i},ligand);
end
