function edges = meshconnee(t) 
%   This function establishes all edges of a triangular mesh given its
%   t-array. The result will be sorted. 
%
%   Copyright (C) 2004-2012 Per-Olof Persson. (See DISTMESH)

    edges      = [t(:,[1, 2]); t(:,[1, 3]); t(:,[2, 3])];         %  All edges duplicated
    edges      = unique(sort(edges, 2), 'rows');                  %  Unique edges as sorted node pairs
end