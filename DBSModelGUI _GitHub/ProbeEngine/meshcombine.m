function [P, t, normals] = meshcombine(varargin)
%   SYNTAX
%   [P, t] = meshcombine(P1, P2, P3, t1, t2, t3, n1, n2, n3);
%   DESCRIPTION
%   This function combines multiple meshes into one mesh 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    MN = size(varargin, 2); %   total number of arguments
    %   Combine arrays of vertices
    P = [];
    for m = 1:MN/3
        P = [P; varargin{m}];
    end   
    %   Combine arrays of triangles
    temp = 0;
    t = [];
    for m = MN/3+1:2*MN/3        
        t = [t; (varargin{m}+temp)];
        temp = temp + size(varargin{m-MN/3}, 1);
    end
    %   Combine arrays of normal vectors
    normals = [];
    for m = 2*MN/3+1:MN
        normals = [normals; varargin{m}];
    end
end

