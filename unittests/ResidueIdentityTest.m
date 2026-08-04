%% Missing segid resolves only through a unique FASTA identity
residue_index = build_residue_index([1 2 3 4], 'AAAA', {'A','A','A','A'});
[chain, segid, resnum] = resolve_residue_identity('A', '', 2, residue_index, 'fixture.txt', 7);
assert(strcmp(chain, 'A'));
assert(strcmp(segid, 'A'));
assert(resnum == 2);

%% Explicit segid must match exactly
residue_index = build_residue_index([1 2 3 4], 'AAAA', {'A','A','A','A'});
[chain, segid, resnum] = resolve_residue_identity('A', 'A', 3, residue_index, 'fixture.txt', 8);
assert(strcmp(chain, 'A'));
assert(strcmp(segid, 'A'));
assert(resnum == 3);

try
    resolve_residue_identity('A', 'B', 3, residue_index, 'fixture.txt', 9);
    assert(false, 'Expected an explicit segid mismatch to fail.');
catch err
    assert(strcmp(err.identifier, 'RiboDraw:ResidueIdentityNotFound'));
end

%% Missing and ambiguous identities fail closed
residue_index = build_residue_index([1 2 3 4], 'AAAA', {'A','A','A','A'});
try
    resolve_residue_identity('A', '', 99, residue_index, 'fixture.txt', 10);
    assert(false, 'Expected an absent residue identity to fail.');
catch err
    assert(strcmp(err.identifier, 'RiboDraw:ResidueIdentityNotFound'));
end

ambiguous_index = build_residue_index([1 1], 'AA', {'X','Y'});
try
    resolve_residue_identity('A', '', 1, ambiguous_index, 'fixture.txt', 11);
    assert(false, 'Expected an ambiguous residue identity to fail.');
catch err
    assert(strcmp(err.identifier, 'RiboDraw:AmbiguousResidueIdentity'));
end

%% Duplicate exact FASTA identities are rejected
try
    build_residue_index([1 1], 'AA', {'X','X'});
    assert(false, 'Expected duplicate exact residue identities to fail.');
catch err
    assert(strcmp(err.identifier, 'RiboDraw:DuplicateResidueIdentity'));
end

%% All annotation readers use the same canonical residue identities
residue_index = build_residue_index([1 2 3 4], 'AAAA', {'A','A','A','A'});
testdata = fullfile(fileparts(mfilename('fullpath')), 'testdata');
stems = read_stems(fullfile(testdata, 'segid_identity_case.stems.txt'), residue_index);
base_pairs = read_base_pairs(fullfile(testdata, 'segid_identity_case.base_pairs.txt'), residue_index);

assert(length(stems) == 1);
assert(all(strcmp(stems{1}.segid1, {'A','A'})));
assert(all(strcmp(stems{1}.segid2, {'A','A'})));
assert(length(base_pairs) == 2);
assert(strcmp(base_pairs{1}.segid1, 'A'));
assert(strcmp(base_pairs{1}.segid2, 'A'));

stem_pair_keys = validate_stem_annotations(stems, base_pairs);
assert(length(stem_pair_keys) == 2);

%% Legacy direct calls without a FASTA index retain their old parse result
testdata = fullfile(fileparts(mfilename('fullpath')), 'testdata');
legacy_stems = read_stems(fullfile(testdata, 'segid_identity_case.stems.txt'));
legacy_pairs = read_base_pairs(fullfile(testdata, 'segid_identity_case.base_pairs.txt'));
assert(all(strcmp(legacy_stems{1}.segid1, {'',''})));
assert(isempty(legacy_pairs{1}.segid1));
