function unitnormals = meshnormals(P, t);
%   SYNTAX
%   unitnormals = meshnormals(P, t)
%   DESCRIPTION
%   For a manifold mesh P, t, this function returns an array of outer
%   normal vectors for each triangle
%   Inputs:
%   P - Vertex array of the mesh (N x 3)
%   t - Facets of the mesh (N x 3)


vert1 = P(t(:, 1),:);
vert2 = P(t(:, 2),:);
vert3 = P(t(:, 3),:);
% Finding edges
edge1 = vert2 - vert1;
edge2 = vert3 - vert1;

normal = cross(edge1, edge2, 2);                                 % Calculating the normal of the triangle
length = sqrt(normal(:,1).^2+normal(:,2).^2+normal(:,3).^2);     % Calculating length of the normal
unitnormals= normal./(repmat(length,size(normal,3),3));          % Normalization of the normal vector
