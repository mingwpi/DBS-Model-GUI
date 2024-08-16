function E = bemf5_volume_field_electric(Points, c, P, t, Center, Area, normals, R, planeABCD)
%   Computes electric field for an array Points anywhere in space (line,
%   surface, volume). This field is due to surface charges at triangular
%   facets only. Includes accurate neighbor triangle integrals for
%   points located close to a charged surface.   
%   R is the dimensionless radius of the precise-integration sphere
%
%   Copyright SNM 2017-2020
%   R = is the local radius of precise integration in terms of average triangle size
    

    if(nargin < 9)
        planeABCD = [];
    end
    
    %   FMM 2019
    tic
    srcinfo.sources = Center';                      %   source points
    targ            = Points';                      %   target points
    prec            = 1e-2;                         %   precision    
    pg      = 0;                                    %   nothing is evaluated at sources
    pgt     = 2;                                    %   field and potential are evaluated at target points
    srcinfo.charges = c.'.*Area';                   %   charges
    U               = lfmm3d(prec, srcinfo, pg, targ, pgt);
    E               = -U.gradtarg'/(4*pi);  
    disp([newline 'Plain FMM time ' num2str(toc) ' s']);

    tic    
    %   Undo the effect of the m-th triangle charge on neighbors and
    %   add precise integration instead  
    %   Contribution of the charge of triangle m to the field at all points is sought        
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
    disp([newline 'Neighbor search time ' num2str(toc) ' s']);
    
    timeFullSubtraction = 0;
    timeFullIntegral = 0;
    for j = 1:length(eligibleTriangles)
        index = ineighborlocal{j};
        m = eligibleTriangles(j);
        if ~isempty(index)
        tic
            temp        = repmat(Center(m, :), length(index), 1) - Points(index, :);   %   these are distances to the observation points
            DIST        = sqrt(dot(temp, temp, 2));                                    %   single column                                            %Fast calculation for distance^2
            I           = Area(m)*temp./repmat(DIST.^3, 1, 3);                         %   center-point integral, standard format
            E(index, :) = E(index, :) - (- c(m)*I/const);
            timeFullSubtraction = timeFullSubtraction + toc;            
            tic
            r1      = P(t(m, 1), :);    %   row
            r2      = P(t(m, 2), :);    %   row
            r3      = P(t(m, 3), :);    %   row           
            I       = potint2(r1, r2, r3, normals(m, :), Points(index, :));     %   analytical precise integration MATLAB
            E(index, :)= E(index, :) + (- c(m)*I/const);
            timeFullIntegral = timeFullIntegral + toc;
        end    
    end
    disp([newline 'Plain FMM subtraction time ' num2str(timeFullSubtraction) ' s']); 
    disp([newline 'Adding accurate integrals time ' num2str(timeFullIntegral) ' s']); 
end