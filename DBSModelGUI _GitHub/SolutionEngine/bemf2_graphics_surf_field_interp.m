function [ ] = bemf2_graphics_surf_field_interp(parent, P, t, FQ, Indicator, objectnumber) 
%   Surface field graphics:  plot a field quantity FQ at the surface of a
%   brain compartment with the number "tissuenumber". Interpolates over
%   triangles
%
%   Copyright SNM 2017-2021

    %%  Interpolation for nodes
    tic
    t0              = t(Indicator==objectnumber, :);   
    [Pobs, tobs]    = fixmesh(P, t0);
    TR              = triangulation(tobs, Pobs);
    V               = vertexAttachments(TR);
    N               = size(Pobs, 1); 
    tempnodes       = zeros(N, 1);
    for m = 1:N
        tempnodes(m) = mean(FQ(V{m}));
    end
  
    %%   Graphics - interpolation plot
    Ptemp = Pobs'; ttemp = tobs';
    X = reshape(Ptemp(1, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Y = reshape(Ptemp(2, ttemp(1:3, :)),[3, size(ttemp, 2)]);
    Z = reshape(Ptemp(3, ttemp(1:3, :)),[3, size(ttemp, 2)]); 
    %%   Interpolate field for vertexes - global  
    temp = tempnodes';
    C = temp(ttemp(1:3, :)); 
    %%  Plot
    patch(parent, X, Y, Z, C, 'FaceAlpha', 1.0, 'EdgeColor', 'none', 'FaceLighting', 'flat'); 
    colorbar(parent);
    % colormap hsv;
    % camlight('headlight');
    % lighting phong;
    axis(parent,'equal');  axis(parent,'tight'); 
    xlabel(parent,'x, mm'); ylabel(parent,'y, mm'); zlabel(parent,'z, mm');
    set(parent,'Color','White'); 
end