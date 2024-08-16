function [Int] = potint2(r1, r2, r3, normal, ObsPoint)
%   This function computes potential integrals grad(1/r) for a single triangle
%   The integrals are not divided by the area
%   Vectorized for an arbitrary number of observation points
%
%   Copyright SNM 2004-2020

%   Test (comparison with Wang et. al., 2003):
% %     clear all;
% %     format long;
% %     r1 = [62.5 25.0 0];
% %     r2 = [62.5 25.0 2];
% %     r3 = [62.5 37.5 0];
% %     ObsPoint = [62.5 0.0 0.0];
% %     tempv           = cross(r2-r1, r3-r1);  %   correct normal sign!
% %     temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
% %     normal          = tempv/temps;
% %     Area            = temps/2;   
% %     %   Analytical integral
% %      [coeff, weights, IndexF] = tri(25, 10);
% %      Int = [0 0 0];
% %      for m = 1:length(coeff)     
% %          Point = coeff(1, m)*r1 + coeff(2, m)*r2 + coeff(3, m)*r3;
% %          R = ObsPoint - Point;
% %          Int = Int - Area*weights(m)*R/(sqrt(sum(R.*R)))^3;
% %      end
% %      Int
% %      Int = potint2(r1, r2, r3, normal, ObsPoint)

    N           = size(ObsPoint, 1);
    I           = zeros(N, 3);
    S           = zeros(N, 9);
    Int         = zeros(N, 3);
    BetaTerm    = zeros(N, 3);

    %   Create r+ and r- coordinates
    temp = [r2 r1 r3 r1 r3 r2];
    r    = repmat(temp, N, 1);

    %   Create l coordinates
    normabsl1 = sqrt(sum((r2-r1).^2));
    normabsl2 = sqrt(sum((r3-r1).^2));
    normabsl3 = sqrt(sum((r3-r2).^2));
    temp = [(r2-r1)./normabsl1 (r3-r1)./normabsl2 (r3-r2)./normabsl3];
    l    = repmat(temp, N, 1);

    %   Create unit normal to the edges of the triangle
    u(1:3) = +cross((r2-r1)./normabsl1, normal);
    u(4:6) = -cross((r3-r1)./normabsl2, normal);
    u(7:9) = +cross((r3-r2)./normabsl3, normal);
    u      = repmat(u, N, 1);

    %   Create projection vector of the observation point
    NORM = repmat(normal, N, 1);
    ndot = sum(ObsPoint.*NORM, 2);
    p    = ObsPoint - repmat(ndot, 1, 3).*NORM;

    %   Is the projection point on the triangle edge or its continuation?
    %   Calculate midpoints of the edges of a triangle
    temp     = 0.5*[r1+r2 r1+r3 r2+r3];
    midpoint = repmat(temp, N, 1);
    %   Calculate vectors from observation point project to midpoints
    vector = midpoint - [p p p];
      %   Move the observation  point (projection) from the edge
    factor = 1e-6;
    factor1 = factor*normabsl1;
    factor2 = factor*normabsl2;
    factor3 = factor*normabsl3;
    check1 = sum(vector(:, 1:3).*u(:, 1:3), 2);
    check2 = sum(vector(:, 4:6).*u(:, 4:6), 2);
    check3 = sum(vector(:, 7:9).*u(:, 7:9), 2);	
    index  = (abs(check1)<factor1); ObsPoint(index, :) =  ObsPoint(index, :) - factor1*u(index, 1:3);    
    index  = (abs(check2)<factor2); ObsPoint(index, :) =  ObsPoint(index, :) - factor2*u(index, 4:6);
    index  = (abs(check3)<factor3); ObsPoint(index, :) =  ObsPoint(index, :) - factor3*u(index, 7:9);
    ndot   = sum(ObsPoint.*NORM, 2);
    p      = ObsPoint - repmat(ndot, 1, 3).*NORM;

    % Calculation of the analytical formula
    count = 0;
    for c1 = 0:2    
        %   Distance of observation point perpendicular to the plane with triangle
        temporary   = ObsPoint - r(:, [1:3] + 3*(count + 1));
         %%  Changes compared to potint.m   
        distanceobs = sum(NORM.*temporary, 2);     
         %%  End of changes compared to potint.m   

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

        %  Changes compared to potint.m
        %   A value of one term of the analytic sum 1/R         
        d1      = atan(P0.*lplus./(R0.*R0 + abs(distanceobs).*RPLUS));
        d2      = atan(P0.*lminus./(R0.*R0 + abs(distanceobs).*RMINUS));
        d3      = log((RPLUS + lplus)./(RMINUS + lminus));

        % End of changes compared to potint.m

        dotPHATu= sum(PHAT.*u(:, [1:3] + 3*c1), 2);        

        %  Changes compared to potint.m        
        BetaTerm(:, c1+1) = dotPHATu.*(d1 - d2);
        S(:, 1 + 3*c1) = d3.*u(:, 1 + 3*c1);
        S(:, 2 + 3*c1) = d3.*u(:, 2 + 3*c1);
        S(:, 3 + 3*c1) = d3.*u(:, 3 + 3*c1);
        % End of changes compared to potint.m
        count = count + 2;      
    end

    %  Changes compared to potint.m         
    Beta = sum(BetaTerm, 2);

    I(:, 1) = S(:, 1) + S(:, 4) + S(:, 7);
    I(:, 2) = S(:, 2) + S(:, 5) + S(:, 8);
    I(:, 3) = S(:, 3) + S(:, 6) + S(:, 9);

    %   find sign of distanceobs
    Sign = sign(distanceobs); 

    %   value of integral for 1/R
    Int(:, 1) = -normal(1)*Sign.*Beta - I(:, 1);
    Int(:, 2) = -normal(2)*Sign.*Beta - I(:, 2);
    Int(:, 3) = -normal(3)*Sign.*Beta - I(:, 3);

    %   Contribution is zero when the projection point on the edge (Wilton et al. 1984, p. 279)    
    Int(isnan(Int)) = 0;  Int(isinf(Int)) = 0;    
    % End of changes compared to potint.m
end