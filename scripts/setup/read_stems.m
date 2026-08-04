function helices = read_stems(helix_file, residue_index)
% READ_STEMS Read RiboDraw stem annotations.
%
% helices = read_stems(file)
% helices = read_stems(file, residue_index)
%
% A FASTA-derived residue index makes chain/segid/resnum resolution strict.

helices = {};
if ~exist(helix_file, 'file')
    return;
end
if nargin == 1
    residue_index = [];
elseif nargin ~= 2
    error('RiboDraw:InvalidReaderArguments', ...
        'read_stems expects 1 or 2 input arguments.');
end
fid = fopen(helix_file);
if fid == -1
    error('RiboDraw:AnnotationFileOpenFailed', ...
        'Could not open stem file %s.', helix_file);
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
    if length(cols) < 2
        error('RiboDraw:InvalidStemLine', ...
            '%s:%d expected at least 2 columns.', helix_file, line_number);
    end

    clear helix
    [helix.resnum1, helix.chain1, helix.segid1, ok1] = get_resnum_from_tag(cols{1});
    [helix.resnum2, helix.chain2, helix.segid2, ok2] = get_resnum_from_tag(cols{2});
    if ~ok1 || ~ok2
        error('RiboDraw:InvalidResidueIdentity', ...
            '%s:%d contains an invalid stem range.', helix_file, line_number);
    end
    if length(helix.resnum1) ~= length(helix.resnum2)
        error('RiboDraw:UnequalStemStrands', ...
            '%s:%d stem strands have unequal lengths (%d and %d).', ...
            helix_file, line_number, length(helix.resnum1), length(helix.resnum2));
    end

    if ~isempty(residue_index)
        for i = 1:length(helix.resnum1)
            [helix.chain1(i), helix.segid1{i}, helix.resnum1(i)] = ...
                resolve_residue_identity(helix.chain1(i), helix.segid1{i}, ...
                helix.resnum1(i), residue_index, helix_file, line_number);
            [helix.chain2(i), helix.segid2{i}, helix.resnum2(i)] = ...
                resolve_residue_identity(helix.chain2(i), helix.segid2{i}, ...
                helix.resnum2(i), residue_index, helix_file, line_number);
        end
    end

    if length(cols) > 2
        helix.name = cols{3};
    else
        helix.name = '';
        warning('RiboDraw:UnnamedStem', ...
            'No stem name found for %s/%s in file %s.', cols{1}, cols{2}, helix_file);
    end
    helix.helix_tag = sanitize_tag(sprintf('Helix_%s%s%d', ...
        helix.chain1(1), helix.segid1{1}, helix.resnum1(1)));
    helices = [helices, helix];
end
