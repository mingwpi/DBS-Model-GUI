function t = meshreorient(P, t, normals) 
%   This function reorient triangles (needs to be improved)

%   Copyright SNM 2020
    N           = size(t, 1);
    for m = 1:N
        Vertexes        = P(t(m, 1:3)', :)';
        r1              = Vertexes(:, 1);
        r2              = Vertexes(:, 2);
        r3              = Vertexes(:, 3);
        tempv           = cross(r2-r1, r3-r1);  %   definition (*)
        temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
        normalcheck     = tempv'/temps;
        if sum(normalcheck.*normals(m, :))<0;   %   rearrange vertices to have exactly the outer normal
            t(m, 2:3) = t(m, 3:-1:2);           %   by definition (*)
        end     
    end   
end