function [color, countXZ, PofXZ, EofXZ, TofXZ, pointsXZ, x, z] = bem03_define_plane(bem02out, bem03in)

    %   This script defines three observation planes and plots cross-sections and
    %   NIfTI data when availble
    %
    %   Copyright SNM/WAW 2017-2020

    % From structure to individual fields
    cellfun(@(f) assignin('caller', f, bem02out.(f)), fieldnames(bem02out));
    cellfun(@(f) assignin('caller', f, bem03in.(f)), fieldnames(bem03in));
    
    %% Create AUX grid
    numPoints = 200;
    
    x = linspace(xmin, xmax, numPoints);
    z = linspace(zmin, zmax, numPoints);
    [xx, zz] = meshgrid(x, z);
    xx = reshape(xx, [], 1);
    zz = reshape(zz, [], 1);
    pointsXZ = [xx zz];
    pointsXZ(:, 3) = pointsXZ(:, 2);
    pointsXZ(:, 2) = Y;
    
    %%   Process cross-section data to enable fast (real time) display 
    %   This block finds all edges and attached triangles for separate brain
    %   compartments. This script is required for subsequent visualizations.
    %   Process surface model data
    tic
    %   Preallocate cell arrays
    m_max = 2; %length(name);
    tS = cell(m_max, 1);
    nS = tS; %  Reuse this empty cell array for other initialization
    eS = tS;
    TriPS = tS;
    TriMS = tS;
    PS = P;
    for m = 1:m_max
        tS{m} = t(Indicator == m, :);
        nS{m} = normals(Indicator == m, :);
        [eS{m}, TriPS{m}, TriMS{m}] = mt(tS{m}); 
    end
    SurfaceDataProcessTime = toc
    
    color   = prism(m_max); color(4, :) = [0 1 1];
    PofXZ = cell(m_max, 1);   %   intersection nodes for a tissue
    EofXZ = cell(m_max, 1);   %   edges formed by intersection nodes for a tissue
    TofXZ = cell(m_max, 1);   %   intersected triangles
    NofXZ = cell(m_max, 1);   %   normal vectors of intersected triangles
    countXZ = [];   %   number of every tissue present in the slice
    for m = 1:m_max 
        [Pi, ti, polymask, flag] = meshplaneintXZ(PS, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Y);
        if flag % intersection found                
            countXZ               = [countXZ m];
            PofXZ{m}            = Pi;               %   intersection nodes
            EofXZ{m}            = polymask;         %   edges formed by intersection nodes
            TofXZ{m}            = ti;               %   intersected triangles
            NofXZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
        end
    end
  
end