function [EC] = meshneighborints_En(P, t, normals, Area, Center, RnumberE, ineighborE)
%   Accurate integration for electric field/electric potential on neighbor facets
%   Copyright WAW 2020
    N = size(t, 1);
    integralxc      = zeros(RnumberE, N);    %   center-point Ex integrals for array of neighbor triangles 
    integralyc      = zeros(RnumberE, N);    %   center-point Ey integrals for array of neighbor triangles 
    integralzc      = zeros(RnumberE, N);    %   center-point Ez integrals for array of neighbor triangles 
    
    gauss       = 25;   %   number of integration points in the Gaussian quadrature  
                        %   for the outer potential integrals
                        %   Numbers 1, 4, 7, 13, 25 are permitted 
    %   Gaussian weights for analytical integration (for the outer integral)
    if gauss == 1;  [coeffS, weightsS, IndexS]  = tri(1, 1); end
    if gauss == 4;  [coeffS, weightsS, IndexS]  = tri(4, 3); end
    if gauss == 7;  [coeffS, weightsS, IndexS]  = tri(7, 5); end
    if gauss == 13; [coeffS, weightsS, IndexS]  = tri(13, 7); end
    if gauss == 25; [coeffS, weightsS, IndexS]  = tri(25, 10); end
    if gauss == 0;  [coeffS, weightsS, IndexS]  = tri(50);      end
    

    %   Main loop for analytical double integrals (parallel, 24 workers)
    %   This is the loop over columns of the system matrix
    tic
    if(isempty(gcp('nocreate')))
        error('A parallel pool must be initialized prior to running meshneighborints_En');
    end
    integrale = zeros(N, RnumberE);
    parfor n = 1:N                  %   inner integral; (n =1 - first column of the system matrix, etc.)        
        % Calculate observation points on this triangle
        ObsPoints = zeros(IndexS, 3);
        for p = 1:IndexS
            ObsPoints(p, :)  = coeffS(1, p)*P(t(n, 1), :) +  coeffS(2, p)*P(t(n, 2), :) +  coeffS(3, p)*P(t(n, 3), :);
        end
        
        % Get vertices of neighbor triangles acting on this triangle
        index = ineighborE(:,n);
        r1 = P(t(index, 1), :); %get first vertex of each neighbor triangle
        r2 = P(t(index, 2), :); %get second vertex of each neighbor triangle
        r3 = P(t(index, 3), :); %get third vertex of each neighbor triangle
        
        %Int_temp stores the contribution of each triangle (column) to each observation point (row)
        Int_temp = potint4b(r1, r2, r3, ObsPoints);
        Int_temp(:,1) = 0; %kill self-term
        %Now weight and sum each column of Int_temp properly to get a single row
        %weightsS: row vector containing contribution of each observation point to final triangle
        Int = weightsS*Int_temp; % Exploiting dimensions of weightsS and Int_temp to ensure proper product occurs
        integrale(n, :) = Int;       
           
        %   Center-point electric-field integrals
        temp    = repmat(Center(n, :), RnumberE, 1) - Center(index, :); %   these are distances to the observation/target triangle
        DIST    = sqrt(dot(temp, temp, 2));                             %   single column                
        I       = Area(n)*temp./repmat(DIST.^3, 1, 3);                  %   center-point integral, standard format    
        I(1, :) = 0;                                                    %   self integrals will give zero
        integralxc(:, n) = -I(:, 1);    %   center-point integrals, entries of non-zero rows of n-th column
        integralyc(:, n) = -I(:, 2);    %   center-point integrals, entries of non-zero rows of n-th column
        integralzc(:, n) = -I(:, 3);    %   center-point integrals, entries of non-zero rows of n-th column        
        
    end
    
    %% Properly weight integrale with the self-triangle area instead of the neighbor-triangle area
    area_neighbor = Area(transpose(ineighborE));
    area_self = repmat(Area, 1, RnumberE);
    integrale = integrale .* area_self ./ area_neighbor;
    
    %%  Define useful sparse matrices EC, PC (for GMRES speed up)    
    const           = 1/(4*pi);  
    integralc       = zeros(RnumberE, N);    %   normal integral component for array of neighbor triangles (center point) - to speed up GMRES
    
    for n = 1:N                  %   inner integral; (n =1 - first column of the system matrix, etc.)             
        index = ineighborE(:, n); %   those are non-zero rows of the system matrix for given n 
                               
        integralc(:, n)  =       +(integralxc(:, n).*normals(index, 1) + ...
                                   integralyc(:, n).*normals(index, 2) + ...
                                   integralzc(:, n).*normals(index, 3));
    end 
    
    ii  = ineighborE;
    jj  = repmat([1:N], RnumberE, 1);
    EC  = sparse(ii, jj, const*(-integralc + transpose(integrale)));               %   almost symmetric
    
    MainLoopParallelTimeE = toc
end