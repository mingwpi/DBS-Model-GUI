function [ ] = bemf1_graphics_electrodes(parent,P, t, IndicatorElectrodes, electrodeVoltages, flag) 
%   Electrode plot (with thick edges)
%
%   Copyright SNM 2017-2022

eV          = electrodeVoltages;
enumber     = length(eV);
VV          = (eV-min(eV))/(max(eV)-min(eV));
index       = min(floor(256*VV)+1, 256);
map         = hsv(512);
map         = map(257:512, :);
strge.Color = map(index, :);

for m = 1:enumber
    if flag == 0    % 3D view
        p = patch(parent,'vertices', P, 'faces', t(IndicatorElectrodes==m, :));
        p.FaceColor = strge.Color(m, :);
        p.EdgeColor = 'k';
        p.LineWidth = 0.5;
    end
    if flag == 1    % XY view
        Q = P; Q(:, 3) = [];
        p = patch(parent,'vertices', Q, 'faces', t(IndicatorElectrodes==m, :));
        p.FaceColor = strge.Color(m, :);
        p.EdgeColor = 'k';
        p.LineWidth = 0.5;
    end
    if flag == 2    % XZ view
        Q = P; Q(:, 2) = [];
        p = patch(parent,'vertices', Q, 'faces', t(IndicatorElectrodes==m, :));
        p.FaceColor = strge.Color(m, :);
        p.EdgeColor = 'k';
        p.LineWidth = 0.5;
    end
    if flag == 3    % YZ view
        Q = P; Q(:, 1) = [];
        p = patch(parent,'vertices', Q, 'faces', t(IndicatorElectrodes==m, :));
        p.FaceColor = strge.Color(m, :);
        p.EdgeColor = 'k';
        p.LineWidth = 0.5;
    end
end
