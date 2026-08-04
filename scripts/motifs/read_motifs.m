function motifs = read_motifs( motif_file, residue_index )
%  motifs = read_motifs( motif_file );
%
% read text file with Rosetta rna_motif output format. Lines are like:
%
% U_TURN A:33-35
% U_TURN A:55-57
% UA_HANDLE A:53-54 A:58 A:61
% T_LOOP A:53-58 A:61
% INTERCALATED_T_LOOP A:53-58 A:61 A:18
%
%  INPUT
%   motifs_file = text file with Rosetta rna_motif output.
%
% (C) R. Das, Stanford University, 2019

motifs = {};
if ~exist( motif_file, 'file' ) return; end;
if nargin == 1
    residue_index = [];
end
fid = fopen( motif_file );
if fid == -1
    error('RiboDraw:AnnotationFileOpenFailed', ...
        'Could not open motif file %s.', motif_file);
end
cleanup = onCleanup(@() fclose(fid));
line_number = 0;
while ~feof( fid )
    line = fgetl( fid );
    line_number = line_number + 1;
    % UA_HANDLE A:53-54 A:58 A:61
    if ~ischar(line); break; end
    if isempty(strtrim(line)); continue; end
    cols = strsplit(strtrim(line));
    if length(cols) < 2 && ~isempty(residue_index)
        error('RiboDraw:InvalidMotifLine', ...
            '%s:%d expected at least 2 columns.', motif_file, line_number);
    end
    if length( cols ) >= 2 
        clear motif

        motif.motif_type = cols{1};
        [resnum,chains,segid] = get_resnum_from_tag( strjoin(cols(2:end)) );
        motif.associated_residues = {};
        for i = 1:length( resnum )
            if ~isempty(residue_index)
                [chains(i),segid{i},resnum(i)] = resolve_residue_identity( ...
                    chains(i),segid{i},resnum(i),residue_index,motif_file,line_number);
            end
            motif.associated_residues = [ motif.associated_residues, sprintf( 'Residue_%s%s%d',  chains(i), segid{i}, resnum(i) ) ];            
        end
        assert( length( motif.associated_residues ) > 1 );
        motif.motif_tag = sprintf( 'Motif_%s_%s', motif.motif_type, strrep(motif.associated_residues{1},'Residue_','') );
        motifs = [motifs,motif];
    end
end


