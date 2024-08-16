function [En] = bemf4_surface_field_electric_accurate(c, Center, Area, normals, EC)    
%   Computes the normal electric field just inside
%
%   Copyright SNM 2017-2020

    %tic
    [P0, E0]      = bemf4_surface_field_electric_plain(c, Center, Area);        %   Plain FMM result    
    correction  = EC*c;                                                         %   Correction of plain FMM result
    En          = +c/2 + correction ...                                         %   This is the dominant (exact) matrix part and the "undo" terms for center-point FMM
                       + sum(normals.*E0, 2);                                   %   This is the full center-point FMM part      
end
