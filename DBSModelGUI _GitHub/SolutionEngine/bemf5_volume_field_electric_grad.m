function [Ex, Ey, Ez] = bemf5_volume_field_electric_grad(Points, c, P, t, Center, Area, normals, R)
%   Computes electric field tensor gradient for an array Points anywhere in space (line,
%   surface, volume). This tensor is due to surface charges at triangular
%   facets only. Includes accurate neighbor triangle integrals for
%   points located close to a charged surface.   
%   R is the dimensionless radius of the precise-integration sphere
%
%   Copyright SNM 2017-20204
%   R = is the local radius of precise integration in terms of average triangle size
    
    %   FMM 2019
    srcinfo.nd      = 3;                            %   three charge vectors    
    srcinfo.sources = Center';                      %   source points
    targ            = Points';                      %   target points
    prec            = 1e-2;                         %   precision    
    pg      = 0;                                    %   nothing is evaluated at sources
    pgt     = 2;                                    %   potential/gradient is evaluated at target points
    
    %   set of dipoles
    const                    = 1/(4*pi);  % normalization
    temp                     = c.'.*Area';
    srcinfo.dipoles(1, 1, :) = 1*temp;
    srcinfo.dipoles(1, 2, :) = 0*temp;
    srcinfo.dipoles(1, 3, :) = 0*temp;
    srcinfo.dipoles(2, 1, :) = 0*temp;
    srcinfo.dipoles(2, 2, :) = 1*temp;
    srcinfo.dipoles(2, 3, :) = 0*temp;
    srcinfo.dipoles(3, 1, :) = 0*temp;
    srcinfo.dipoles(3, 2, :) = 0*temp;
    srcinfo.dipoles(3, 3, :) = 1*temp;   
    U                        = lfmm3d(prec, srcinfo, pg, targ, pgt);
    Ex                       = const*squeeze(U.gradtarg(1, :, :)).';
    Ey                       = const*squeeze(U.gradtarg(2, :, :)).';
    Ez                       = const*squeeze(U.gradtarg(3, :, :)).';
    
    %   Undo the effect of the m-th triangle charge on neighbors and
    %   add zero instead  
    %   Contribution of the charge of triangle m to the field at all points is sought        
    M = size(Center, 1);      
    const = 4*pi;    
    Size  = mean(sqrt(Area));
    ineighborlocal   = rangesearch(Points, Center, R*Size, 'NSMethod', 'kdtree'); % over triangles: M by X  
    for m =1:M
        index       = ineighborlocal{m};  % index into points that are close to triangle m   
        if ~isempty(index)
            temp        = repmat(Center(m, :), length(index), 1) - Points(index, :);   %   these are distances to the observation points
            DIST        = sqrt(dot(temp, temp, 2));                                    %   single column     
            I0          = Area(m)./DIST.^3;                                            %   center-point integral, standard format 
            %--------------------------------------------------------------------------------------------------------------------
            Ix          = 3*Area(m)*temp(:, 1).*temp(:, 1)./DIST.^5;                   %   center-point integral, standard format 
            Iy          = 3*Area(m)*temp(:, 1).*temp(:, 2)./DIST.^5;                   %   center-point integral, standard format 
            Iz          = 3*Area(m)*temp(:, 1).*temp(:, 3)./DIST.^5;                   %   center-point integral, standard format             
            Ex(index, 1) = Ex(index, 1) + (- c(m)*I0/const);
            Ex(index, 1) = Ex(index, 1) - (- c(m)*Ix/const);
            Ex(index, 2) = Ex(index, 2) - (- c(m)*Iy/const);
            Ex(index, 3) = Ex(index, 3) - (- c(m)*Iz/const); 
            %--------------------------------------------------------------------------------------------------------------------
            Ix          = 3*Area(m)*temp(:, 2).*temp(:, 1)./DIST.^5;                   %   center-point integral, standard format 
            Iy          = 3*Area(m)*temp(:, 2).*temp(:, 2)./DIST.^5;                   %   center-point integral, standard format 
            Iz          = 3*Area(m)*temp(:, 2).*temp(:, 3)./DIST.^5;                   %   center-point integral, standard format             
            Ey(index, 2) = Ey(index, 2) + (- c(m)*I0/const);
            Ey(index, 1) = Ey(index, 1) - (- c(m)*Ix/const);
            Ey(index, 2) = Ey(index, 2) - (- c(m)*Iy/const);
            Ey(index, 3) = Ey(index, 3) - (- c(m)*Iz/const); 
            %--------------------------------------------------------------------------------------------------------------------
            Ix          = 3*Area(m)*temp(:, 3).*temp(:, 1)./DIST.^5;                   %   center-point integral, standard format 
            Iy          = 3*Area(m)*temp(:, 3).*temp(:, 2)./DIST.^5;                   %   center-point integral, standard format 
            Iz          = 3*Area(m)*temp(:, 3).*temp(:, 3)./DIST.^5;                   %   center-point integral, standard format             
            Ez(index, 3) = Ez(index, 3) + (- c(m)*I0/const);
            Ez(index, 1) = Ez(index, 1) - (- c(m)*Ix/const);
            Ez(index, 2) = Ez(index, 2) - (- c(m)*Iy/const);
            Ez(index, 3) = Ez(index, 3) - (- c(m)*Iz/const); 
        end    
    end     
end