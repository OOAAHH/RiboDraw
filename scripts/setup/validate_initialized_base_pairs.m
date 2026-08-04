function validate_initialized_base_pairs(stems, base_pairs)
% VALIDATE_INITIALIZED_BASE_PAIRS Audit linker classification after setup.

stem_pair_keys = validate_stem_annotations(stems, base_pairs);
stem_pair_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:length(stem_pair_keys)
    stem_pair_set(stem_pair_keys{i}) = true;
end

expected_stem_count = length(stem_pair_keys);
expected_noncanonical_count = length(base_pairs) - expected_stem_count;

for n = 1:length(stems)
    stem = stems{n};
    stem_length = length(stem.resnum1);
    for k = 1:stem_length
        partner_position = stem_length - k + 1;
        residue_tag1 = sanitize_tag(sprintf('Residue_%s%s%d', ...
            stem.chain1(k), stem.segid1{k}, stem.resnum1(k)));
        residue_tag2 = sanitize_tag(sprintf('Residue_%s%s%d', ...
            stem.chain2(partner_position), stem.segid2{partner_position}, ...
            stem.resnum2(partner_position)));
        residue1 = getappdata(gca, residue_tag1);
        residue2 = getappdata(gca, residue_tag2);
        if ~isfield(residue1, 'stem_partner') || ...
                ~strcmp(residue1.stem_partner, residue_tag2) || ...
                ~isfield(residue2, 'stem_partner') || ...
                ~strcmp(residue2.stem_partner, residue_tag1)
            error('RiboDraw:NonreciprocalStemPartner', ...
                'Stem partner relation is not reciprocal for %s / %s.', ...
                residue_tag1, residue_tag2);
        end
    end
end

for i = 1:length(base_pairs)
    base_pair = base_pairs{i};
    residue_tag1 = sanitize_tag(sprintf('Residue_%s%s%d', ...
        base_pair.chain1, base_pair.segid1, base_pair.resnum1));
    residue_tag2 = sanitize_tag(sprintf('Residue_%s%s%d', ...
        base_pair.chain2, base_pair.segid2, base_pair.resnum2));
    pair_key = make_pair_key( ...
        make_residue_key(base_pair.chain1, base_pair.segid1, base_pair.resnum1), ...
        make_residue_key(base_pair.chain2, base_pair.segid2, base_pair.resnum2));

    if isKey(stem_pair_set, pair_key)
        expected_type = 'stem_pair';
    else
        expected_type = 'noncanonical_pair';
    end
    linker_tag = sanitize_tag(sprintf('Linker_%s%s%d_%s%s%d_%s', ...
        base_pair.chain1, base_pair.segid1, base_pair.resnum1, ...
        base_pair.chain2, base_pair.segid2, base_pair.resnum2, expected_type));
    if ~isappdata(gca, linker_tag)
        error('RiboDraw:BasePairClassificationMismatch', ...
            'Missing %s linker for %s / %s.', expected_type, residue_tag1, residue_tag2);
    end

    linker = getappdata(gca, linker_tag);
    if ~strcmp(linker.type, expected_type)
        error('RiboDraw:BasePairClassificationMismatch', ...
            'Pair %s / %s was classified as %s instead of %s.', ...
            residue_tag1, residue_tag2, linker.type, expected_type);
    end
end

actual_stem_count = length(get_tags('Linker_', 'stem_pair'));
actual_noncanonical_count = length(get_tags('Linker_', 'noncanonical_pair'));
if actual_stem_count ~= expected_stem_count || ...
        actual_noncanonical_count ~= expected_noncanonical_count
    error('RiboDraw:BasePairClassificationCountMismatch', ...
        ['Expected %d stem and %d noncanonical pair linkers, ' ...
         'but initialized %d and %d.'], ...
        expected_stem_count, expected_noncanonical_count, ...
        actual_stem_count, actual_noncanonical_count);
end


function key = make_residue_key(chain, segid, resnum)
key = [chain, char(31), segid, char(31), sprintf('%d', resnum)];


function key = make_pair_key(key1, key2)
ordered = sort({key1, key2});
key = [ordered{1}, char(30), ordered{2}];
