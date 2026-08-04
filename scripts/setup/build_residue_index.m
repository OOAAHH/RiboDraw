function residue_index = build_residue_index(resnum, chains, segid)
% BUILD_RESIDUE_INDEX Build the authoritative residue identity index from FASTA.

assert(length(resnum) == length(chains), ...
    'RiboDraw:InvalidResidueIndex', 'resnum and chains must have equal lengths.');
assert(length(resnum) == length(segid), ...
    'RiboDraw:InvalidResidueIndex', 'resnum and segid must have equal lengths.');

residue_index.resnum = reshape(resnum, 1, []);
residue_index.chains = reshape(chains, 1, []);
residue_index.segid = reshape(segid, 1, []);
residue_index.exact_to_position = containers.Map('KeyType', 'char', 'ValueType', 'double');
residue_index.chain_res_to_positions = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:length(residue_index.resnum)
    chain = residue_index.chains(i);
    residue_number = residue_index.resnum(i);
    segment_id = residue_index.segid{i};
    exact_key = make_exact_key(chain, segment_id, residue_number);
    display_tag = make_display_tag(chain, segment_id, residue_number);

    if isKey(residue_index.exact_to_position, exact_key)
        error('RiboDraw:DuplicateResidueIdentity', ...
            'FASTA contains duplicate residue identity %s.', display_tag);
    end
    residue_index.exact_to_position(exact_key) = i;

    chain_res_key = make_chain_res_key(chain, residue_number);
    if isKey(residue_index.chain_res_to_positions, chain_res_key)
        positions = residue_index.chain_res_to_positions(chain_res_key);
        residue_index.chain_res_to_positions(chain_res_key) = [positions, i];
    else
        residue_index.chain_res_to_positions(chain_res_key) = i;
    end
end


function key = make_exact_key(chain, segid, resnum)
key = [chain, char(31), segid, char(31), sprintf('%d', resnum)];


function key = make_chain_res_key(chain, resnum)
key = [chain, char(31), sprintf('%d', resnum)];


function tag = make_display_tag(chain, segid, resnum)
if isempty(segid)
    tag = sprintf('%s:%d', chain, resnum);
else
    tag = sprintf('%s:%s:%d', chain, segid, resnum);
end
