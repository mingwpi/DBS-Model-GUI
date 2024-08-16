function P = meshlaplace(P, t, nodes, alpha)
%   SYNTAX 
%   P = meshlaplace3D(P, t, nodes, alpha)
%   DESCRIPTION 
%   This function implements lumped Laplacian smoothing 
%   based on existing Delaunay connectivity for a given set of nodes. Use
%   nodes = 1:size(P, 1) for the global smoothing
%   Inputs:
%   P - array of vertices; t - array of faces; nodels - nodes to be moved;
%   alpha - weighting parameter
%  
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    %   Triangles attached to every vertex (cell array)
    [si] = meshconnvt(t);
    Pnew = P;
    for m = nodes                
        index = unique(reshape(t(si{m}, :), 1, 3*length(si{m})));
        index(index==m) = [];           %   Indexes into vertices connected to vertex m               
        Pnew(m, :) = alpha*sum(P(index, :))/length(index)+(1-alpha)*P(m, :);   
    end
    P = Pnew;
end

function [si] = meshconnvt(t)
%   SYNTAX 
%   [si] = meshconnvt(t)
%   DESCRIPTION 
%   This function returns triangles attached to every vertex (a cell array si) 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

   N  = max(max(t));
    si = cell(N, 1);
    for m = 1:N
        temp  = (t(:,1)==m)|(t(:,2)==m)|(t(:,3)==m);
        si{m} = find(temp>0);
    end 
end

