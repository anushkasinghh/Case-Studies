% FOUND ON THE PAGE https://www.mathworks.com/matlabcentral/answers/1567423-how-to-write-multiple-object-surfaces-solids-to-a-single-stl-file#answer_1567955

function [F V] = t2fv(T,V0)
% FV = T2FV(T)
% [F V] = T2FV(T)
% [F V] = T2FV(FV)
% FV = T2FV(F,V) 
% Simple function to grab F,V data from a triangulation or trirep object.
% I made this because i hate dealing with verbose and inconsistent naming.
% I just want to swap easily and concisely between data formats.
%
% Input can be any of the following:
%    triangulation or delaunayTriangulation objects
%    TriRep or DelaunayTri objects
%    FV structs
%    or separate F,V lists.
% 
% Output can be independent F,V lists or an FV struct.

switch nargin
	case 1
		switch lower(class(T))
			case {'triangulation','delaunaytriangulation'}
				F = T.ConnectivityList;
				V = T.Points;
			case {'trirep','delaunaytri'}
				F = T.Triangulation;
				V = T.X;
			case 'struct'
				% these names are consistent with surf2patch(), reducepatch(), and similar.
				% i'm just going to assume the field names are common.
				% surely there's no way that'll ever blow up in my face. 
				F = T.faces;
				V = T.vertices;
			otherwise
				error('T2FV: expected either a triangulation or trirep object, or an FV struct')
		end
		if nargout == 1
			F = struct('faces',F,'vertices',V);
		end
		
	case 2
		% it's handy but i hate to document something that ruins the name
		F = struct('faces',T,'vertices',V0);
		
	otherwise
		error('T2FV: incorrect number of arguments')
end






