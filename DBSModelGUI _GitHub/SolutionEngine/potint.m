function [I, IRho] = potint(r1, r2, r3, normal, ObsPoint)
%   This function computes single potential integrals I = Is(1/r) and 
%   integrals IRho=Is(vec(r)/r), respectively,  and
%   follows the method described in: Wilton DR, Rao SM, Glisson AW,
%   Schaubert DH, Al-Bundak OM, Butler CM. Potential integrals for uniform
%   and linear source distribution on polygonal and polyhedral domains.
%   IEEE Trans Antennas Propag 1984;AP-32 (3):276–281 - uses the last
%   formula in Eq. (5) on page 279.
%   Inputs:
%   r1, r2, r3 - three triangle vertices (1x3 rows)
%   normal - unit normal vector to the triangle:
%             tempv           = cross(r2-r1, r3-r1);
%             temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
%             normal          = tempv'/temps; 
%   ObsPoint - Nx3 array of observation points
%   Outputs:
%   I - Nx1 array of integrals I = Is(1/r)
%   IRho - Nx3 array of integrals IRho = Is(vec(r)/r)
%   To obtain the complete integrals these results should be multiplied by 
%   triangle area A
%
%   Copyright SNM 2002-2020

    N = size(ObsPoint, 1);
    
    Z = zeros(N, 9);
    S = zeros(N, 3);

    %   Create r+ and r- coordinates
    temp = [r2 r1 r3 r1 r3 r2];
    r    = repmat(temp, N, 1);

    %   Create l coordinates
    normabsl1 = sqrt(sum((r2-r1).^2));
    normabsl2 = sqrt(sum((r3-r1).^2));
    normabsl3 = sqrt(sum((r3-r2).^2));
    temp0 = [(r2-r1)./normabsl1 (r3-r1)./normabsl2 (r3-r2)./normabsl3];
    l    = repmat(temp0, N, 1);

    %   Create unit normal to the edges of the triangle
    u(1:3) = +cross((r2-r1)./normabsl1, normal);
    u(4:6) = -cross((r3-r1)./normabsl2, normal);
    u(7:9) = +cross((r3-r2)./normabsl3, normal);
    u      = repmat(u, N, 1);

    %   Create projection vector of the observation point
    NORM = repmat(normal, N, 1);
    ndot = sum(ObsPoint.*NORM, 2);
    p    = ObsPoint - repmat(ndot, 1, 3).*NORM;
    
    % Calculation of the analytic formula
    count = 0;
    for c1 = 0:2    
        %   Distance of observation point perpendicular to the plane with triangle
        temporary   = ObsPoint - r(:, [1:3] + 3*(count + 1));
        distanceobs = abs(sum(NORM.*temporary, 2));     

        %   Calculate p+ and p-
        d1 = sum(NORM.*r(:, [1:3] + 3*count), 2);
        d2 = sum(NORM.*r(:, [1:3] + 3*(count + 1)), 2);

        pplus  = r(:, [1:3] + 3*count)       - NORM.*repmat(d1, 1, 3);
        pminus = r(:, [1:3] + 3*(count + 1)) - NORM.*repmat(d2, 1, 3);

        %   Calculate l+ and l-
        lplus  = sum(l(:, [1:3] + 3*c1).*(pplus - p), 2);
        lminus = sum(l(:, [1:3] + 3*c1).*(pminus - p), 2);

        %   Perpendicular distance from projection vector to edge
        P0 = abs( sum(u(:, [1:3] + 3*c1).*(pminus-p), 2) ); 

        %   Distances to l+ and l- from projection vector
        PPLUS  = sqrt(P0.*P0 + lplus.*lplus);
        PMINUS = sqrt(P0.*P0 + lminus.*lminus);

        %   Vector containing line P0 measured
        PHAT = (pminus - p - repmat(lminus, 1, 3).*l(:, [1:3] + 3*c1))./repmat(P0, 1, 3);       

        %   Distances to l+ and l- from observation point
        RPLUS  = sqrt(PPLUS.*PPLUS + distanceobs.*distanceobs);
        RMINUS = sqrt(PMINUS.*PMINUS + distanceobs.*distanceobs);
        R0	   = sqrt(P0.*P0 + distanceobs.*distanceobs);

        %   A value of one term of the analytic sum 1/R
        d1      = P0.*log((RPLUS+lplus)./(RMINUS+lminus));        
        d2      = distanceobs.*(atan((P0.*lplus)./(R0.*R0 + distanceobs.*RPLUS)));        
        d3      = distanceobs.*(atan((P0.*lminus)./(R0.*R0 + distanceobs.*RMINUS)));
        
        dotPHATu= sum(PHAT.*u(:, [1:3] + 3*c1), 2);
        S(:, c1+1)   = (d1 - d2 + d3).*dotPHATu;

        %   A value of one term of the analytic sum RHO /R
        d1 = R0.*R0.*log((RPLUS + lplus)./(RMINUS + lminus));
        d2 = lplus.*RPLUS - lminus.*RMINUS;
        d3 = d1 + d2;
        Z(:, [1:3] + c1*3) = repmat(d3, 1, 3).*u(:, [1:3] + c1*3);
        
        count = count + 2;
    end
    
    %   Contribution is zero when the projection point on the edge (Wilton et al. 1984, p. 279)
    S(isnan(S)) = 0;  S(isinf(S)) = 0;  
    Z(isnan(Z)) = 0;  Z(isinf(Z)) = 0; 
    I       = sum(S, 2);    
    IRho    = 0.5*(Z(:, 1:3) + Z(:, 4:6) + Z(:, 7:9));   
end
