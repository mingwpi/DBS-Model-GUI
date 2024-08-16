function p = bemf5_volume_field_potential(Points, c, P, t, Center, Area, normals, R, planeABCD)
%   Computes electric potential for an array Points anywhere in space (line,
%   surface, volume). This field is due to surface charges at triangular
%   facets only. Includes accurate neighbor triangle integrals for
%   points located close to a charged surface.   
%   R is the dimensionles radius of the precise-integration sphere
%
%   Copyright SNM 2018-2020
%   R = is the local radius of precise integration in terms of average triangle size
    
    if(nargin < 9)
        planeABCD = [];
    end

    %   FMM 2019
    srcinfo.sources = Center';                      %   source points
    targ            = Points';                      %   target points
    prec            = 1e-2;                         %   precision    
    pg      = 0;                                    %   nothing is evaluated at sources
    pgt     = 1;                                    %   potential is evaluated at target points
    srcinfo.charges = c.'.*Area';                   %   charges
    U               = lfmm3d(prec, srcinfo, pg, targ, pgt);
    p               = +U.pottarg'/(4*pi);  
    
    if flag == 0    % Only the center-point approximation is used
        return;
    end
    %   Contribution of the charge of triangle m to the field at all points is sought                 
    %   Undo the effect of the m-th triangle charge on neighbor obs. points and
    %   add precise integration instead 
    M = size(Center, 1); 
    const = 4*pi;     
    Size  = mean(sqrt(Area));
    if(isempty(planeABCD))
        eligibleTriangles = 1:size(t, 1);
    else
        d1 = abs(planeABCD(1)*Center(:,1) + planeABCD(2)*Center(:,2) + planeABCD(3)*Center(:,3) + planeABCD(4));
        d2 = norm(planeABCD(1:3));
        d = d1./d2;
        eligibleTriangles = find(d <= R*Size);
    end
    ineighborlocal   = rangesearch(Points, Center(eligibleTriangles, :), R*Size, 'NSMethod', 'kdtree'); % over triangles: M by X  
    %for m =1:M
    for j = 1:length(eligibleTriangles)
        index = ineighborlocal{j};
        m = eligibleTriangles(j);
        if ~isempty(index)
            temp        = repmat(Center(m, :), length(index), 1) - Points(index, :);   %   these are distances to the observation points
            DIST        = sqrt(dot(temp, temp, 2));                                    %   single column                
            I           = Area(m)./DIST;                                               %   center-point integral, standard format    
            p(index)    = p(index) - c(m)*I/const;         
            r1      = P(t(m, 1), :);    %   row
            r2      = P(t(m, 2), :);    %   row
            r3      = P(t(m, 3), :);    %   row           
            I       = potint(r1, r2, r3, normals(m, :), Points(index, :));     %   analytical precise integration MATLAB            
            p(index)= p(index) + c(m)*I/const;
        end    
    end                        
end

