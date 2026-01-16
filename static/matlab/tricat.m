% FOUND ON THE PAGE https://www.mathworks.com/matlabcentral/answers/1567423-how-to-write-multiple-object-surfaces-solids-to-a-single-stl-file#answer_1567955

function [F V] = tricat(Tc)
%  [F V] = TRICAT(TC)
%  T = TRICAT(TC)
%  Combine multiple sets of triangulated surface data into a single representation.
%
%  TC is a cell array, the elements of which may be any of the following:
%    triangulation or delaunayTriangulation objects
%    TriRep or DelaunayTri objects
%    FV structs
%    2-element cell arrays, each containing a set of F,V lists.
%
%  Output is either a single set of F,V lists or a triangulation object.
%    If the user is running R2012b or older, TriRep will be used instead
%    of a triangulation object.
%
%  When combining delaunay-type objects, the constraints are not kept.

	% i originally went to a whole heap of effort to avoid the array growing
	% but this is actually consistently faster for me.  it's certainly much simpler.
	F = []; V = [];
	for k = 1:numel(Tc)
		switch lower(class(Tc{k}))
			case {'triangulation','trirep','delaunaytriangulation','delaunaytri','struct'}
				[Fk Vk] = t2fv(Tc{k});
			case 'cell'
				Fk = Tc{k}{1};
				Vk = Tc{k}{2};
			otherwise
				error(['TRICAT: TC{%d} is not an accepted type. ' ...
					'Elements of TC are expected to be either ' ...
					'triangulation, delaunayTriangulation, or TriRep objects ' ...
					'or either FV structs or 2-element cell arrays containing F and V lists.'],k)
		end
		F = [F; Fk+size(V,1)]; %#ok<*AGROW>
		V = [V; Vk];
	end
	
	% prune the vertex list
	[F V] = prunedupeverts(F,V);

	% if only 1 output is requested, it's a triangulation object.
	% if we're using an old version without triangulation(), use TriRep()
	if nargout == 1
		if verLessThan('matlab','8.1')
			F = TriRep(F,V); %#ok<DTRIREP>
		else
			F = triangulation(F,V);
		end
	end
end