function LHS = bemf4_surface_field_lhs_v(c, Center, Area, Contrast, normals, M, EC, PC, ElectrodeIndexes, electrodeVoltages, weight)   
%   Computes the left hand side of the charge equation for surface charges
%   for active, floating, and sheet electrodes
%
%   Copyright SNM 2023-2024
%   This algorithm may not be used for commercial purposes 

    %   LHS is the user-defined function of c equal to c - Z_times_c which is
    %   exactly the left-hand side of the matrix equation Zc = b
    [P0, E0]      = bemf4_surface_field_electric_plain(c, Center, Area);    %   Plain FMM result    
    correction  = EC*c;                                                     %   Correction of plain FMM result
    LHS         = +c - 2*correction ...                                     %   This is the dominant (exact) matrix part and the "undo" terms for center-point FMM
                     - 2*(Contrast.*sum(normals.*E0, 2));                   %   This is the full center-point FMM part
    
    %   This is the modification for voltage electrodes
    correctionP  = PC*c;                                                    %   Correction of plain FMM result for potential
    P            = P0 + correctionP;                                        %   Exact results for potential

    %   Normal field just OUTSIDE the electrode, both for active, floating,
    %   and sheet electrodes
    En          = +c/2 + correction ...                                     %   This is the dominant (exact) matrix part and the "undo" terms for center-point FMM
                       + sum(normals.*E0, 2);                               %   This is the full center-point FMM part     

    %   Loops over all electrodes
    indexfull   = [];
    indexactive = [];
    PCORR       = [];

    %   Compute PCORR, indexactive
    for j = 1:length(electrodeVoltages)                                      %   LHS for potential with local preconditioners, electrode by electrode
        indexv      = ElectrodeIndexes{j};
        indexfull   = [indexfull; indexv];
        if isnan(electrodeVoltages(j))
            temp        = P(indexv) - mean(P(indexv));
            PCORR       = [PCORR; temp];
        else
            indexactive = [indexactive; indexv];
            PCORR       = [PCORR; P(indexv)];
        end
    end

    %   Apply preconditioner to all electrodes at once
    LHS(indexfull) = M*PCORR;

    % %   Apply current correction to all active electrodes at once
    % I = sum(En(indexactive).*Area(indexactive))/sum(Area(indexactive));
    % LHS(indexactive) = LHS(indexactive) + weight*I;                          %   zero current for all active electrodes combined

    %   Apply current correction to floating electrodes individually
    for j = 1:length(electrodeVoltages)                                      %   LHS for potential with local preconditioners, electrode by electrode
        indexv      = ElectrodeIndexes{j};
        if isnan(electrodeVoltages(j))
            I           = sum(En(indexv).*Area(indexv))/sum(Area(indexv)); 
            LHS(indexv) = LHS(indexv) + weight*I;                           %   zero current for every floating electrode individually
        end
    end
                                                                           
end
