function C = meshtricenter(P, t) 
%   SYNTAX
%   C = meshtricenter(P, t) 
%   DESCRIPTION
%   This function returns triangle centers (a Nx3 array)
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    C   = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));
end
   
    
