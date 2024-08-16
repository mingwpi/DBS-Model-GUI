function bem02out = bem02_mesh_combiner(bem01out)

    %   This is a mesh processor script: it computes basis triangle parameters
    %   and necessary potential integrals, and constructs a combined mesh of a
    %   multi-object structure (for example, a head or a whole body)
    
    %%  Define EM constants
    eps0        = 8.85418782e-012/1e3;  %   Dielectric permittivity of vacuum(~air) F/mm
    mu0         = 1.25663706e-006/1e3;  %   Magnetic permeability of vacuum(~air) H/mm
    
    %% Identify tissue filenames
    clear name;
    name{1} = 'geometry.mat';           %   Load volume first
    load(name{1});
    name{2} = 'probe.mat';              %   Load the probe next
    
    %%  Load tissue meshes and combine individual meshes into a single mesh
    PP = [];
    tt = [];
    nnormals = [];
    Indicator = [];
    IndicatorElectrodesTemp = IndicatorElectrodes;
    %   Combine individual meshes into a single mesh
    tic
    for m = 1:length(name)
        if m == 1
            load(name{m});
        else
            % From structure to individual fields
            cellfun(@(f) assignin('caller', f, bem01out.(f)), fieldnames(bem01out));
        end
        t = meshreorient(P, t, normals);
        tt = [tt; t+size(PP, 1)];
        PP = [PP; P];
        nnormals                = [nnormals; normals];    
        Indicator               = [Indicator; repmat(m, size(t, 1), 1)];
        if m > 1
            IndicatorElectrodes(IndicatorElectrodes>0) = IndicatorElectrodes(IndicatorElectrodes>0) + max(IndicatorElectrodesTemp);
            IndicatorElectrodesTemp = [IndicatorElectrodesTemp; IndicatorElectrodes];
        end
        disp(['Successfully loaded file [' name{m} ']']);
    end
    t = tt;
    P = PP;
    normals = nnormals;
    IndicatorElectrodes = IndicatorElectrodesTemp;
    disp([newline 'Base meshes loaded in ' num2str(toc) ' s']);
    
    %% Put electrodes up front sequentially (1, 2, 3, etc.)
      %  Now put electrodes up front sequentially (1, 2, 3, etc.)
    tt = [];
    nn = [];
    ie = [];
    ii = [];
    for m = 1:max(IndicatorElectrodes) 
        index  = IndicatorElectrodes==m;
        tt     = [tt; t(index, :)];
        nn     = [nn; normals(index, :)];
        ie     = [ie; m*ones(sum(index), 1)];
        ii     = [ii; Indicator(index)];
    end
    index               = IndicatorElectrodes==0;
    t                   = [tt; t(index, :)];
    normals             = [nn; normals(index, :)];
    Indicator           = [ii; Indicator(index)];
    IndicatorElectrodes = [ie; IndicatorElectrodes(index)];
    
    %%  Assign conductivities
    M           = size(t, 1);
    condWM      = 1e-4;     % S/mm 
    
    condin      = zeros(M, 1);
    condout     = zeros(M, 1);
    Contrast    = zeros(M, 1);
    
    %   Lead
    condout(Indicator==1)   = condWM;
    condin(Indicator==1)    = 0;
    condout(Indicator==2)   = condWM;
    condin(Indicator==2)    = 0;
    
    %  Lead, at electrodes (this is critical)
    condin(IndicatorElectrodes>0)   = 1e12;
    
    %   Contrast
    Contrast = (condin-condout)./(condin+condout);
    
    %%  Fix triangle orientation (just in case, optional)
    tic
    t = meshreorient(P, t, normals);
    
    %%   Process other mesh data
    Center      = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));  %   face centers
    Area        = meshareas(P, t);  
    disp([newline 'Triangle properties computed in ' num2str(toc) ' s']);
    
    %%  Check for and process triangles that have coincident centroids
    tic
    disp('Checking combined mesh for duplicate facets ...');
    [P, t, normals, Center, Area, Indicator, condin, condout, Contrast] = ...
        clean_coincident_facets(P, t, normals, Center, Area, Indicator, condin, condout, Contrast);
    disp('Resolved all duplicate facets');
    N           = size(t, 1);
    disp([newline 'Duplicate facets resolved in ' num2str(toc) ' s']);

    % From individual fields to structure
    bem02out = struct('strge', strge, 'P', P, 't', t, 'normals', normals, 'Indicator', Indicator, 'Center', Center, 'Area', Area, ...
                      'condin', condin, 'condout', condout, 'Contrast', Contrast, 'IndicatorElectrodes', IndicatorElectrodes,...
                      'eps0', eps0, 'mu0', mu0);

end