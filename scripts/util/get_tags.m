function tags = get_tags( headstring, tailstring, objnames )
%  tags = get_tags( headstring )
%  tags = get_tags( headstring, tailstring )
%  tags = get_tags( headstring, tailstring, objnames )
%
% Find all objects in figure with tags like:
%
%  headstring...
%
%  or if tailstring is specified:
%
%  headstring...tailstring
%
% INPUTS
%
%   headstring = string that should residue at beginning of tag, e.g., 'Residue_'
%   tailstring = string that should residue at end of tag, e.g., 'domain'. 
%                  [Default: '']
%   objnames   = list of tags to look through [Default: all objects in current drawing
%                  (gca)]
%
% (C) R. Das, Stanford University, 2017

if ~exist( 'objnames','var' )
    vals = getappdata( gca );
    objnames = fields( vals );
end
if ~exist( 'tailstring', 'var' ); tailstring = ''; end
tags = {};
tag_ok = zeros(length(objnames),1);
for n = 1:length( objnames )

    % 只接受以 headstring 开头的 tag，避免误把例如
    %  'TertiaryContact_Residue_A23_Residue_A26' 之类对象当作 Residue_
    %  相关的 tag。
    if length( objnames{n} ) < length( headstring ); continue; end;
    if ~strcmp( objnames{n}(1:length(headstring)), headstring ); continue; end;

    if ~isempty( tailstring )
        if length( objnames{n} ) < length(tailstring ) | ...
                ~strcmp( objnames{n}(end-length(tailstring)+1 : end), tailstring ); continue; end;
    end
    tag_ok(n) = 1;
end
tags = objnames( find( tag_ok ) );
