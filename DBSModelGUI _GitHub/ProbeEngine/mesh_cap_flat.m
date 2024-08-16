function [Pcap, tcap, ncap] = mesh_cap_flat(startend, Pcenter, P, t, M)

    if startend == 2
        PathVector      = Pcenter(end, :) - Pcenter(end-1, :);
        Pend    = P(end-M+1:end, :);
        offset  = Pcenter(end, :);
        Pendrel = Pend - repmat(offset, M, 1);
    end
    if startend == 1
        PathVector      = Pcenter(1, :) - Pcenter(2, :);
        Pend    = P(1:M, :);
        offset  = Pcenter(1, :);
        Pendrel = Pend - repmat(offset, M, 1);
    end
    UnitPathVector  = PathVector/norm(PathVector);
    %   Fill out the planar terminal surface with the nodes
    par     = linspace(0.1, 0.90, round(M/4));
    rad     = par;
    Pcap    = Pendrel;
    for m = 1:length(rad)
        Pcap = [Pcap; rad(m)*Pendrel(1:round(1/rad(m)):end, :)];
    end
    CapNondes = size(Pcap, 1);
    %   Do triangulation and create the (planar) mesh
    X  = max(Pcap(:, 1)) - min(Pcap(:, 1));
    Y  = max(Pcap(:, 2)) - min(Pcap(:, 2));
    Z  = max(Pcap(:, 3)) - min(Pcap(:, 3));
    if Z<=X&Z<=Y
        dt  = delaunayTriangulation(Pcap(:, [1 2]));  %   2D Delaunay
    end
    if Y<=X&Y<=Z
        dt  = delaunayTriangulation(Pcap(:, [1 3]));  %   2D Delaunay
    end
    if X<=Z&X<=Y
        dt  = delaunayTriangulation(Pcap(:, [2 3]));  %   2D Delaunay
    end
    tcap = dt.ConnectivityList;
    %   Perform Laplacian smoothing
    alpha = 0.5;
    nodes = M+1:size(Pcap, 1);
    Pcap                = meshlaplace(Pcap, tcap, nodes, alpha);
    %   Add offset
    Pcap = Pcap + repmat(offset, size(Pcap, 1), 1);
    %   Find the normal vectors
    ncap                = meshnormals(Pcap, tcap);
    index               = dot(ncap, repmat(PathVector, size(tcap, 1), 1), 2)<0;
    ncap(index, :)      = - ncap(index, :);
    