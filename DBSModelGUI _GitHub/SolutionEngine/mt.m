function [edges, TriP, TriM] = mt(t)
%   This function finds all edges and attached triangles of a triangular
%   mesh. It is fast.
%   Michael Tsuk, MathWorks, 2109
    edgesdup   = [t(:,[1, 2]); t(:,[1, 3]); t(:,[2, 3])];         %  All edges duplicated
    e1         = sort(edgesdup, 2);                               %  All edges now in same order 
    [edges, eia(:,1)] = unique(e1, 'rows');                       %  eia(:,1) contains the row in e1 that was chosen for edges
    [~,eia(:,2)] = unique(e1, 'last', 'rows');                    %  eia(:,2) contains the other row in e1 for the same edges
    e2t        = mod(eia-1,size(t,1))+1;                          %  since "edgesdup" contains three copies of t, mod makes the translation.
    e2t        = sort(e2t,2);                                     %  Put the lowest index triangle into TriP.
    TriP       = e2t(:,1);
    TriM       = e2t(:,2);
end

