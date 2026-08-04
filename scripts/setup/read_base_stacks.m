function base_stacks = read_base_stacks( base_stacks_file, residue_index )
% base_stacks = read_base_stacks( base_stacks_file )
%
%  Read .base_stacks.txt file output by Rosetta rna_motif executable
%      should include all base-base combinations in which one resides in the
%      appropriate 'z' range above/below the other one.
%
% INPUT
%
%  base_stacks_file = text file with lines like
%
%                      A:1  B:20 A P
%
%                    i.e.,
%
%                      chain1[:segid1]:resnum1 chain2[:segid2]:resnum2  side orientation
%
%                    where side        is A/B (for the second nucleobase being above/below first nucleobase)
%                          orientation is A/P (antiparallel/parallel based on how the two base's normal vectors
%                                                   are aligned)
%
%
% OUTPUT
%
%  base_stacks       = cell of struct()s with the same information. Reordered so that 
%                         the residue that has an earlier chain/segid (or if same, earlier resnum) is first.
%
%
% See also: SETUP_BASE_STACKS, READ_BASE_PAIRS
% 
% (C) R. Das, Stanford University, 2017

base_stacks = {};
if ~exist( base_stacks_file, 'file' ) return; end;
if nargin == 1
    residue_index = [];
end
fid = fopen( base_stacks_file );
if fid == -1
    error('RiboDraw:AnnotationFileOpenFailed', ...
        'Could not open stack file %s.', base_stacks_file);
end
cleanup = onCleanup(@() fclose(fid));
line_number = 0;
while ~feof( fid )
    line = fgetl( fid );
    line_number = line_number + 1;
    if ~ischar(line) || isempty(strtrim(line)); continue; end
    % C:1347 C:1599 W W C 
    cols = strsplit(strtrim(line));
    if length(cols) < 4
        if isempty(residue_index); continue; end
        error('RiboDraw:InvalidStackLine', ...
            '%s:%d expected 4 columns.', base_stacks_file, line_number);
    end
        [base_stack.resnum1,base_stack.chain1,base_stack.segid1,ok1] = get_one_resnum_from_tag( cols{1} );
        [base_stack.resnum2,base_stack.chain2,base_stack.segid2,ok2] = get_one_resnum_from_tag( cols{2} );
        if ~ok1 || ~ok2
            error('RiboDraw:InvalidResidueIdentity', ...
                '%s:%d contains an invalid residue tag.', base_stacks_file, line_number);
        end
        if ~isempty(residue_index)
            [base_stack.chain1,base_stack.segid1,base_stack.resnum1] = resolve_residue_identity( ...
                base_stack.chain1,base_stack.segid1,base_stack.resnum1,residue_index,base_stacks_file,line_number);
            [base_stack.chain2,base_stack.segid2,base_stack.resnum2] = resolve_residue_identity( ...
                base_stack.chain2,base_stack.segid2,base_stack.resnum2,residue_index,base_stacks_file,line_number);
        end
        base_stack.side = cols{3};
        base_stack.orientation = cols{4};
        base_stacks = [base_stacks, base_stack];
end
