% FOUND ON THE PAGE https://www.mathworks.com/matlabcentral/answers/1567423-how-to-write-multiple-object-surfaces-solids-to-a-single-stl-file#answer_1567955

function [F,V] = prunedupeverts(F,V,tol)
% [F,V] = PRUNEDUPEVERTS(F,V,{TOL})
% Remove duplicate vertices from FV data, remapping F as needed
% 
% See also: prunebadverts, pruneunusedverts

% remove duplicate vertices
if nargin<3
	[V,~,ic] = unique(V,'rows');
	F = ic(F);
else
	[V,~,ic] = uniquetol(V,tol,'byrows',true);
	F = ic(F);
end


