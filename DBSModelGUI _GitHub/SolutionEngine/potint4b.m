% This function calculates n*grad(1/r) at a given observation point obsPoint
% given a triangle with vertices r1, r2, and r3 and normal vector normal.
% It uses the solid-angle approximation of Van Oosterom and Strackee 1983
% to quickly compute the normal component of the field in the vicinity of
% the triangle.
% r1: Nx3 first triangle vertex location for N triangles
% r2: Nx3 second triangle vertex location for N triangles
% r3: Nx3 third triangle vertex location for N triangles
% obsPoint: Mx3 list of observation points at which electric field should be evaluated
% Int: MxN matrix of integral contributions to each point.  Right-multiply by column
%  vector of triangle weights (e.g. charges) to obtain total contribution to each
%  observation point.

% Copyright William Wartman 2020

function [Int] = potint4b(r1, r2, r3, obsPoint)
    %Vectorize operation for triangles and observation points simultaneously
     N = size(r1, 1);        %N triangles
     M = size(obsPoint, 1);  %M observation points
     
    %Dimension 1: triangle index.  Dimension 2: 3. Dimension 3: Observation point index
    r1Exp = repmat(r1, 1, 1, M);
    r2Exp = repmat(r2, 1, 1, M);
    r3Exp = repmat(r3, 1, 1, M);

    obsPointExpA = zeros(1, 3, M);
    obsPointExpA(1, :, :) = transpose(obsPoint);
    obsPointExp = repmat(obsPointExpA, N, 1, 1);

    %Vectors from observation points to triangle vertices (N by 3 by M)
    R1 = r1Exp - obsPointExp;
    R2 = r2Exp - obsPointExp;
    R3 = r3Exp - obsPointExp;
    
    %Norms of vectors (N by 1 by M)
    R1norm = vecnorm(R1, 2, 2);
    R2norm = vecnorm(R2, 2, 2);
    R3norm = vecnorm(R3, 2, 2);
    
    %N by 1 by M
    numerator = dot(R1, cross(R2, R3, 2), 2);
    
    %N by 1 by M
    denominator = (R1norm .* R2norm .* R3norm) + R3norm .* dot(R1, R2, 2) + R2norm .* dot(R1, R3, 2) + R1norm .* dot(R2, R3, 2);
    
    omega = 2*atan2(numerator, denominator);
    
    %Squeeze out the middle dimension and shape omega properly to be scaled
    %by a column vector of triangle charges
    if(N ~= 1) 
        Int = transpose(squeeze(omega));
    else
        % If there is only one triangle, the first dimension is 1 and thus
        % is squeezed out along with the second dimension, leaving a row 
        % vector as desired.
        Int = squeeze(omega);
    end
       
end
