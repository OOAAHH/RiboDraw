function [chain, segid, resnum] = resolve_residue_identity(chain, segid, resnum, residue_index, source_file, line_number)
% RESOLVE_RESIDUE_IDENTITY Resolve one annotation residue against FASTA.

if ~exist('source_file', 'var') || isempty(source_file)
    source_file = '<unknown>';
end
if ~exist('line_number', 'var') || isempty(line_number)
    line_number = 0;
end

assert(isscalar(chain) && isscalar(resnum), ...
    'RiboDraw:InvalidResidueIdentity', ...
    '%s:%d has an invalid residue identity.', source_file, line_number);

if ~isempty(segid)
    exact_key = make_exact_key(chain, segid, resnum);
    if ~isKey(residue_index.exact_to_position, exact_key)
        error('RiboDraw:ResidueIdentityNotFound', ...
            '%s:%d residue %s is not present in FASTA.', ...
            source_file, line_number, make_display_tag(chain, segid, resnum));
    end
    position = residue_index.exact_to_position(exact_key);
else
    chain_res_key = make_chain_res_key(chain, resnum);
    if ~isKey(residue_index.chain_res_to_positions, chain_res_key)
        error('RiboDraw:ResidueIdentityNotFound', ...
            '%s:%d residue %s is not present in FASTA.', ...
            source_file, line_number, make_display_tag(chain, segid, resnum));
    end
    positions = residue_index.chain_res_to_positions(chain_res_key);
    if length(positions) ~= 1
        candidates = cell(1, length(positions));
        for i = 1:length(positions)
            candidates{i} = make_display_tag( ...
                residue_index.chains(positions(i)), ...
                residue_index.segid{positions(i)}, ...
                residue_index.resnum(positions(i)));
        end
        error('RiboDraw:AmbiguousResidueIdentity', ...
            '%s:%d residue %s has multiple FASTA matches: %s.', ...
            source_file, line_number, make_display_tag(chain, segid, resnum), ...
            strjoin(candidates, ', '));
    end
    position = positions(1);
end

chain = residue_index.chains(position);
segid = residue_index.segid{position};
resnum = residue_index.resnum(position);


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
