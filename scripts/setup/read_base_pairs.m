function base_pairs = read_base_pairs(base_pairs_file, residue_index, chains_ref, segid_ref)
% READ_BASE_PAIRS Read RiboDraw base-pair annotations.
%
% base_pairs = read_base_pairs(file)
% base_pairs = read_base_pairs(file, residue_index)
% base_pairs = read_base_pairs(file, resnum, chains, segid)
%
% When a FASTA-derived residue index is supplied, every residue identity is
% resolved strictly against it. Missing segid values are completed only for
% a unique (chain,resnum) match.

base_pairs = {};
if ~exist(base_pairs_file, 'file')
    return;
end

if nargin == 1
    residue_index = [];
elseif nargin == 4
    residue_index = build_residue_index(residue_index, chains_ref, segid_ref);
elseif nargin ~= 2
    error('RiboDraw:InvalidReaderArguments', ...
        'read_base_pairs expects 1, 2, or 4 input arguments.');
end

fid = fopen(base_pairs_file);
if fid == -1
    error('RiboDraw:AnnotationFileOpenFailed', ...
        'Could not open base-pair file %s.', base_pairs_file);
end
cleanup = onCleanup(@() fclose(fid));

line_number = 0;
while ~feof(fid)
    line = fgetl(fid);
    line_number = line_number + 1;
    if ~ischar(line) || isempty(strtrim(line))
        continue;
    end
    cols = strsplit(strtrim(line));
    if length(cols) < 5
        error('RiboDraw:InvalidBasePairLine', ...
            '%s:%d expected 5 columns.', base_pairs_file, line_number);
    end

    [base_pair.resnum1, base_pair.chain1, base_pair.segid1, ok1] = ...
        get_one_resnum_from_tag(cols{1});
    [base_pair.resnum2, base_pair.chain2, base_pair.segid2, ok2] = ...
        get_one_resnum_from_tag(cols{2});
    if ~ok1 || ~ok2
        error('RiboDraw:InvalidResidueIdentity', ...
            '%s:%d contains an invalid residue tag.', base_pairs_file, line_number);
    end

    if ~isempty(residue_index)
        [base_pair.chain1, base_pair.segid1, base_pair.resnum1] = ...
            resolve_residue_identity(base_pair.chain1, base_pair.segid1, ...
            base_pair.resnum1, residue_index, base_pairs_file, line_number);
        [base_pair.chain2, base_pair.segid2, base_pair.resnum2] = ...
            resolve_residue_identity(base_pair.chain2, base_pair.segid2, ...
            base_pair.resnum2, residue_index, base_pairs_file, line_number);
    end

    base_pair.edge1 = cols{3};
    base_pair.edge2 = cols{4};
    base_pair.LW_orientation = cols{5};
    base_pairs = [base_pairs, ordered_base_pair(base_pair)];
end

