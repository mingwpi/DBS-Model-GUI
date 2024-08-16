function bemf2_graphics_vol_field(parent,temp, th1, th2, levels, a, b)
%   Volume field graphics:  plot a field quantity temp in the observation
%   plane using a planar contour plot. Revision 071318
%
%   temp - quantity to plot
%   th1, th2 - two threshold levels introduced manually
%   levels - number of levels in the contour plot introduced manually
%   a, b - x and y arguments
%
%   Copyright SNM 2017-2020

    temp(temp>+th1) = +th1;
    temp(temp<+th2) = +th2;     
    [C, h]          = contourf(parent,a, b, reshape(temp, length(a), length(b)), levels);
    tick            = round((th1-th2)/levels, 1, 'significant');
    h.LevelList     = tick*round(h.LevelList/tick);
    h.ShowText      = 'on';
    h.LabelSpacing  = 50;
    clabel(C, h,'FontSize',21,'Color','white')
    colorbar(parent);
end