function bem01out = bem01_DBS_probe_creator(bem01in, strge)
    %   This script creates triangular surface mesh for a long (bent) object
    %   with imprinted electrodes
    %
    %   Copyright SNM 2017-2024 

    % From structure to individual fields
    cellfun(@(f) assignin('caller', f, bem01in.(f)), fieldnames(bem01in));

    %%  Set paths to local engine folders
    s = pwd;
    if(~isunix) slash = '\'; else slash = '/'; end
    warning off; rmpath(genpath(s)); warning on;
    engine_path =   [s, slash, 'ProbeEngine']; 
    addpath(engine_path);
    engine_path =   [s, slash, 'SolutionEngine']; 
    addpath(engine_path);
    
    %%  Create probe's centerline
    N               = 160;                  %   Subdivisions
    par             = linspace(0, L, N);
    x               = 0*par;                %   Segments
    y               = 0*par;                %   Segments
    z               = par;                  %   Segments
    Pcenter(:, 1)   = x';
    Pcenter(:, 2)   = y';
    Pcenter(:, 3)   = z';
    
    %%  Create object's cross-section: a cylindrical cross section
    R           = a/2;
    M           = 49;           %   number of cross-section subdivisions (better be odd)
    x           = a/2*cos(2*pi*[0:M-1]/M);
    y           = a/2*sin(2*pi*[0:M-1]/M);
    
    %%  Create 2-manifold surface mesh in the form of the (bent) cylinder
    [P, t, normals] = meshsurface(Pcenter, x, y); 
    %%  Create cap(s)
    startend = 1;
    [Pe, te, ne] = mesh_cap_sphere(startend, Pcenter, P, t, M, a, a);
    %[Pe, te, ne] = mesh_cap_flat(startend, Pcenter, P, t, length(x));
    startend = 2;
    [Ps, ts, ns] = mesh_cap_flat(startend, Pcenter, P, t, length(x));
    %%  Combine meshes with caps
    [P, t, normals] = meshcombine(P, Pe, Ps, t, te, ts, normals, ne, ns);
    size(P)
    [P, t] = fixmesh(P, t, 1e-6);
    
    %%  Input data for imprinting electrodes (from GUI)
  
    
    strge.R    = R;                                                 %   already defined above, mm
    
    %%  Imprinting electrodes
    %   Angles are running from 0 to pi and from -pi to 0
    strge.NumberOfElectrodes    = length(strge.type); 
    angle                       = strge.angle;      %   in degrees
    angle2                      = 120-angle;        %   angles must run from -180 to 180 (due to atan)
    strge.phi1 = [];
    strge.phi2 = [];
    k = 0;
    for m = 1:strge.NumberOfElectrodes
        if strge.type(m)==1
            strge.phi1 = [strge.phi1, 0];
            strge.phi2 = [strge.phi2, 0];
        else
            k = k + 1;
            if k == 1       
                strge.phi1 = [strge.phi1, -0.5*angle];
                strge.phi2 = [strge.phi2, +0.5*angle];
            end
            if k == 2
                strge.phi1 = [strge.phi1, 0.5*angle+1*angle2];
                strge.phi2 = [strge.phi2, 1.5*angle+1*angle2];
            end
            if k == 3
                strge.phi1 = [strge.phi1, 1.5*angle+2*angle2-360];
                strge.phi2 = [strge.phi2, 2.5*angle+2*angle2-360];
                k = 0;
            end
        end
    end
    strge.phi1 = strge.phi1/180*pi;
    strge.phi2 = strge.phi2/180*pi;
    
    [P, t, normals, IndicatorElectrodes] = meshimprint_sector(P, t, normals, strge);

    % From individual fields to structure
    bem01out = struct('strge', strge, 'P', P, 't', t, 'normals', normals, 'IndicatorElectrodes', IndicatorElectrodes);
end
 

