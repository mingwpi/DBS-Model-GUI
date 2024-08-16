function bem12out = bem12_charge_engine_it(bem02out, bem11out, electrodeVoltages, ElectrodeIndexes, indexv, V)
    %   This script computes the induced surface charge density for an
    %   inhomogeneous multi-tissue object given the primary electric field, with
    %   accurate neighbor integration
    %
    %   Copyright SNM/WAW 2017-2024
    
    % From structure to individual fields
    cellfun(@(f) assignin('caller', f, bem02out.(f)), fieldnames(bem02out));
    cellfun(@(f) assignin('caller', f, bem11out.(f)), fieldnames(bem11out));
    
    %%  Parameters of the iterative solution
    iter        = 50;                      %    Maximum possible number of inner iterations in the solution 
    ITER        = 4;                       %    Maximum possible number of outer iterations in the solution 
    relres      = 1e-7;                    %    Minimum acceptable relative residual 
    
    %%  Solution for voltage electrodes
    weight          = 1.0;                  %   Global current conservation law 
                                            %   To improve the current
                                            %   conservation accuracy, increase
                                            %   weight (convergence slows down)
    b           = zeros(size(t, 1), 1);     %   Right-hand side of the matrix equation
    b(indexv)   = M*V;                      %   All electrodes held at constant voltage (active, floating, sheet)
    
    %  GMRES iterative solution     
    MATVEC      = @(cv) bemf4_surface_field_lhs_v(cv, Center, Area, Contrast, normals, M, EC, PC, ElectrodeIndexes, electrodeVoltages, weight); 
    tic;
    h   = waitbar(0.5, 'Please wait - Running GMRES');
    [cv, its, resvec] = fgmres(MATVEC, b, relres, 'restart', iter, 'max_iters', ITER, 'x0', b);
    %[cv, flag, rres, its, resvec] = gmres(MATVEC, b, [], relres, iter, [], [], b); 
    GMREStime = toc, close(h);
    %   Total electrode currents
    En = bemf4_surface_field_electric_accurate_En(cv, Center, Area, normals, EC); %    just outside
    electrodeCurrents = zeros(length(ElectrodeIndexes), 1);  
    electrodeImpedances = zeros(length(ElectrodeIndexes), 1);
    for j = 1:length(ElectrodeIndexes)
        index = ElectrodeIndexes{j};
        electrodeCurrents(j) = +sum(condout(index).*En(index).*Area(index));
        electrodeImpedances(j) = electrodeCurrents(j);
    end
    electrodeCurrents

    RESVEC = [];
    for m = 1:size(resvec, 2)
        if m == size(resvec, 2)
            RESVEC = [RESVEC; resvec(1:its(2), m)];
        else
            RESVEC = [RESVEC; resvec(:, m)];
        end
    end

    bem12out = struct('cv', cv, 'electrodeCurrents', electrodeCurrents, 'RESVEC', RESVEC);
end
