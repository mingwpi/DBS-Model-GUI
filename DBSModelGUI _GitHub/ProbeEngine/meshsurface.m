function [P, t, normals] = meshsurface(Pcenter, x, y)
%   Outputs a 2-manifold P, t mesh for a single arbitrarily 
%   bent conductor in 3D. The conductor could be either open or closed. In the
%   last case, the start point and the end point must coincide.
%   Inputs:
%   Pcenter - centerline of the conductor in 3D Pcenter(:, 1:3) in meters;
%   x - x-coordinates of the cross-section contour in the xy-plane in meters;
%   y - y-coordinates of the cross-section contour in the xy-plane in meters;
%   Outputs:
%   P  -  P-aray of surface mesh vertices
%   t  -  t-array of surface triangular facets
%   Copyright SNM 2022
    
    %   Create vector segments along the path
    PathVector      = Pcenter(2:end, :) - Pcenter(1:end-1, :);
    %   Add termination segments for the directional plane method
    Closed          = norm(Pcenter(1, :)-Pcenter(end, :))<1024*eps;
    if Closed
        PathVector = [PathVector; PathVector(1, :)];
    else
        PathVector = [PathVector; PathVector(end, :)];
    end
    %   Create unit directional vectors
    UnitPathVector      = PathVector./repmat(vecnorm(PathVector')', 1, 3);
    %   Align the object cross-section with the first point of the first segment
    %   and perpendicular to the starting path direction (define array ptemp)
    Pcross(:, 1) = x';
    Pcross(:, 2) = y';
    Pcross(:, 3) = 0;
    edges(:, 1) = [1:size(Pcross, 1)  ]';
    edges(:, 2) = [2:size(Pcross, 1) 1]';
    NE          = size(edges, 1);   % number of nodes/edges in the cross-section
    Pcross      = meshrotate1(Pcross, ...
        UnitPathVector(1, 1), UnitPathVector(1, 2), UnitPathVector(1, 3));
    ptemp       = Pcross + repmat(Pcenter(1, :), NE, 1); 
    %   Add nodes/triangles for the side surface (the main loop)
    t     = [];
    P     = ptemp;
    Steps = size(PathVector, 1) - 1;
    for m = 1:Steps
        %   Define the normal vector of the directional plane
        PlaneNormal = UnitPathVector(m, :) + UnitPathVector(m+1, :);
        PlaneNormal = PlaneNormal/norm(PlaneNormal);
        %   Construct intersections with the directional plane
        for n = 1:NE
            [ptemp(n, :),  ~] = line_plane_intersection(UnitPathVector(m, :), ptemp(n, :), PlaneNormal, Pcenter(m+1, :));
        end
        P = [P; ptemp];
        %   Local connectivity: from the previous layer to the next layer
        t1(:, 1:2)      = edges;                        %   Lower nodes        
        t1(:, 3)        = edges(:, 1) + NE;             %   Upper node
        t2(:, 2:-1:1)   = edges       + NE;             %   Upper nodes
        t2(:, 3)        = edges(:, 2);                  %   Lower node
        ttemp           = [t1; t2];
        t               = [t; ttemp+NE*(m-1)];
    end
    normals = meshnormals(P, t);
end

