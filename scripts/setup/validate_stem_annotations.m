function stem_pair_keys = validate_stem_annotations(stems, base_pairs)
% VALIDATE_STEM_ANNOTATIONS Validate stem topology before drawing mutation.

base_pair_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:length(base_pairs)
    base_pair = base_pairs{i};
    key1 = make_residue_key(base_pair.chain1, base_pair.segid1, base_pair.resnum1);
    key2 = make_residue_key(base_pair.chain2, base_pair.segid2, base_pair.resnum2);
    base_pair_set(make_pair_key(key1, key2)) = true;
end

partner_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
stem_pair_keys = {};
for n = 1:length(stems)
    stem = stems{n};
    strand1_length = length(stem.resnum1);
    strand2_length = length(stem.resnum2);
    if strand1_length ~= strand2_length
        error('RiboDraw:UnequalStemStrands', ...
            'Stem %d has unequal strand lengths (%d and %d).', ...
            n, strand1_length, strand2_length);
    end

    for k = 1:strand1_length
        partner_position = strand1_length - k + 1;
        key1 = make_residue_key(stem.chain1(k), stem.segid1{k}, stem.resnum1(k));
        key2 = make_residue_key(stem.chain2(partner_position), ...
            stem.segid2{partner_position}, stem.resnum2(partner_position));

        if strcmp(key1, key2)
            error('RiboDraw:SelfPairedStemResidue', ...
                'Stem %d pairs residue %s with itself.', n, display_key(key1));
        end
        assert_new_partner(partner_map, key1, key2, n);
        assert_new_partner(partner_map, key2, key1, n);

        pair_key = make_pair_key(key1, key2);
        if ~isempty(base_pairs) && ~isKey(base_pair_set, pair_key)
            error('RiboDraw:StemPairMissingFromBasePairs', ...
                'Stem %d pair %s / %s is absent from base_pairs.', ...
                n, display_key(key1), display_key(key2));
        end
        stem_pair_keys{end + 1} = pair_key;
    end
end


function assert_new_partner(partner_map, residue_key, partner_key, stem_number)
if isKey(partner_map, residue_key)
    error('RiboDraw:DuplicateStemPartner', ...
        'Stem %d assigns residue %s more than one stem partner.', ...
        stem_number, display_key(residue_key));
end
partner_map(residue_key) = partner_key;


function key = make_residue_key(chain, segid, resnum)
key = [chain, char(31), segid, char(31), sprintf('%d', resnum)];


function key = make_pair_key(key1, key2)
ordered = sort({key1, key2});
key = [ordered{1}, char(30), ordered{2}];


function tag = display_key(key)
parts = strsplit(key, char(31));
if isempty(parts{2})
    tag = sprintf('%s:%s', parts{1}, parts{3});
else
    tag = sprintf('%s:%s:%s', parts{1}, parts{2}, parts{3});
end
