function [ElectrodeIndexes, indexv, V] = bem04_configure_electrodes(electrodeVoltages, bem02out) 
    %   This script introduces electrode data
    %   Copyright SNM/WAW 2018-2022

    % From structure to individual fields
    cellfun(@(f) assignin('caller', f, bem02out.(f)), fieldnames(bem02out));
    
    %%  Define global electrode indexes (cell array ElectrodeIndexes)
    ElectrodeIndexes = cell(max(IndicatorElectrodes), 1);
    for j = 1:max(IndicatorElectrodes)
        ElectrodeIndexes{j} = find(IndicatorElectrodes==j);
    end
    indexv = transpose(vertcat(ElectrodeIndexes{:}));
    
    %%  Define the voltage excitation vector (floating electrodes must be assigned NaN)
    %   Voltage (V) applied to each electrode facet
    V      = [];   
    for j = 1:max(IndicatorElectrodes)
        if ~isnan(electrodeVoltages(j))
            V = [V; electrodeVoltages(j)*ones(length(ElectrodeIndexes{j}), 1)];
        else
            V = [V; 0.0*ones(length(ElectrodeIndexes{j}), 1)];
        end
    end
end

