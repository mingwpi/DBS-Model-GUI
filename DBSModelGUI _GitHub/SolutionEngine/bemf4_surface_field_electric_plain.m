function [P, E] = bemf4_surface_field_electric_plain(c, Center, Area)
%   Computes potential/continuous electric field on a surface facet due to
%   charges on ALL OTHER facets using plain FMM
%   Self-terms causing discontinuity may not be included
%   To obtain the true field, use E = E/eps0;
%
%   Copyright SNM 2018-2020    
   
    %  FMM 2019   
    %   Fields plus potentials of surface charges (potential not used)
    %   Only FMM (without correction)
    prec            = 1e-4; %   precision-> OK for surfaces
    pg              = 2;    %   potential and field are evaluated
    srcinfo.sources = Center';      %   source/target points
    srcinfo.charges = (c.*Area)';   %   real charges
    U               = lfmm3d(prec, srcinfo, pg);    %   FMM
    P               = +U.pot'/(4*pi);   %   potential
    E               = -U.grad'/(4*pi);  %   field          
end
