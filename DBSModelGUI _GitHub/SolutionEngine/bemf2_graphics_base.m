function [] = bemf2_graphics_base(parent,P, t, c)
%   Surface plot

%   Copyright SNM 2017-2020

    p = patch(parent,'vertices', P, 'faces', t);
    p.FaceColor = c.FaceColor;
    p.EdgeColor = c.EdgeColor;
    p.FaceAlpha = c.FaceAlpha;
    daspect(parent,[1 1 1]);      
    xlabel(parent,'x'); ylabel(parent,'y'); zlabel(parent,'z');
	
    NumberOfTrianglesInShell = size(t, 1);
    edges = meshconnee(t);
    temp = P(edges(:, 1), :) - P(edges(:, 2), :);
    AvgEdgeLengthInShell = mean(sqrt(dot(temp, temp, 2)));
end