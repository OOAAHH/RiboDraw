function other_contacts = read_other_contacts( other_contacts_file, residue_index )
% other_contacts = read_other_contacts( other_contacts_file )
%
%  Read .other_contacts.txt file output by Rosetta rna_motif executable.
%
%    All pairs of nucleotides that make an atom-atom contact less than 3 Angstroms, after
%       filtering out doublets that are recognized as base pairs and base stacks.
%    Mostly hydrogen bonds involving O2' (2' hydroxyl)  and O1P/O2P (phosphate).
%
% INPUT
%
%  other_contacts_file = text file with lines like
%
%                      A:1  B:20 O2' O2'
%
%                    i.e.,
%
%                      chain1[:segid1]:resnum1 chain2[:segid2]:resnum2  atom1 atom2
%
%                    where atom1 and atom2 denote names  of atoms that come within 3 Angstroms.
%
% OUTPUT
%
%  other_contacts       = cell of struct()s with the same information. 
%
%
% See also: READ_BASE_STACKS, READ_BASE_PAIRS.
% 
% (C) R. Das, Stanford University, 2017

other_contacts = {};
if ~exist( other_contacts_file, 'file' ) return; end;
if nargin == 1
    residue_index = [];
end
fid = fopen( other_contacts_file );
if fid == -1
    error('RiboDraw:AnnotationFileOpenFailed', ...
        'Could not open other-contact file %s.', other_contacts_file);
end
cleanup = onCleanup(@() fclose(fid));

line_number = 0;
while ~feof( fid )
    line = fgetl( fid );
    line_number = line_number + 1;
    if ~ischar(line) || isempty(strtrim(line)); continue; end
    % C:1347 C:1599 O2' N3 
    cols = strsplit(strtrim(line));
    if length(cols) < 4
        if isempty(residue_index); continue; end
        error('RiboDraw:InvalidOtherContactLine', ...
            '%s:%d expected 4 columns.', other_contacts_file, line_number);
    end
        [other_contact.resnum1,other_contact.chain1,other_contact.segid1,ok1] = get_one_resnum_from_tag( cols{1} );
        [other_contact.resnum2,other_contact.chain2,other_contact.segid2,ok2] = get_one_resnum_from_tag( cols{2} );
        if ~ok1 || ~ok2
            error('RiboDraw:InvalidResidueIdentity', ...
                '%s:%d contains an invalid residue tag.', other_contacts_file, line_number);
        end
        if ~isempty(residue_index)
            [other_contact.chain1,other_contact.segid1,other_contact.resnum1] = resolve_residue_identity( ...
                other_contact.chain1,other_contact.segid1,other_contact.resnum1,residue_index,other_contacts_file,line_number);
            [other_contact.chain2,other_contact.segid2,other_contact.resnum2] = resolve_residue_identity( ...
                other_contact.chain2,other_contact.segid2,other_contact.resnum2,residue_index,other_contacts_file,line_number);
        end
        other_contact.atom1 = cols{3};
        other_contact.atom2 = cols{4};
        other_contacts = [other_contacts, other_contact];
end
