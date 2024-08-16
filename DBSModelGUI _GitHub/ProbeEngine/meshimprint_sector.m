function [P, t, normals, IndicatorElectrodes] = meshimprint_sector(P, t, normals, strge)

%%   Imprint electrodes
%   Imprints an arbitrary number of electrodes for a DBS probe
%   Electrode surface is an intersection of the host surface with a rectangular cuboid
%   Returns new arrays P, t, normals
%   Returns indexes in t into particular electrodes
%   strge.type = 1 for ring electrodes
%   strge.type = 2 for sectorial electrodes
%   Copyright SNM 2017-2024

    NumberOfElectrodes      = strge.NumberOfElectrodes; 
    
    %%   Establish connectivity
    %   si - triangles attached to every vertex (neighbor triangles)
    %   vi - vertices attached to every vertex (neighbor edges)      
    [si, vi]        = vertices(P, t);
    %   edges - array of mesh edges
    edges           = meshconnee(t);
    %   avgedgelength - average edge length   
    temp            = P(edges(:, 1), :) - P(edges(:, 2), :);                              
    
    %%  Move nodes located close to the boundary exactly to the boundary
    avgedgelength   = mean(sqrt(dot(temp, temp, 2)));
    %   tol - relative tolerance with regard to edge lengt
    tol = 0.5;
    for m = 1:NumberOfElectrodes        
        [DIST, NEARESTPOINT]    = meshimprint_distance(P, strge, m);
        index                   = abs(DIST)/avgedgelength<tol;
        P(index, :)             = NEARESTPOINT(index, :);
    end
    

    %%   Split all intersected triangles    
    %   Add boundary nodes/add triangles for all edges crossing the boundary
    Numbers  = size(P, 1);     %    global numbers of new boundary nodes (accumulating)
    NewTriangles = [];         %    new triangles
    NewNormals   = [];         %    new normal vectors
    NewNodes     = [];         %    new nodes (with repetition)    
    Remove       = [];         %    triangles to remove (to split)    
    for m = 1:NumberOfElectrodes
        %   Find all nodes within the volume
        [DIST, ~]    = meshimprint_distance(P, strge, m);
        temp  = DIST<-1024*eps;
        DIST  = abs(DIST);
        temp  = find(temp>0);   %   all nodes within the boundary    
        for n = 1:length(temp)  %   node, which is in 
            point = temp(n);    %   check this inner node
            index_in        = intersect(vi{point}, temp);   %   all neighbor nodes of the inner node which are inside the sphere
            index_out       = setdiff(vi{point}, temp);     %   all neighbor nodes of the inner node which are outside the sphere                
            if ~isempty(index_out)                          %   find all crossing triangles   
                for p = 1:length(si{point})                 %   loop over triangles attached to the inner node which may intersect the boundary
                    TriNum = si{point}(p);                  %   individual triangle TriNum attached to the inner node which may intersect the boundary
                    IndexIn     = intersect(t(TriNum, 1:3)', index_in);     %    index into node(s) of TriNum that are inside                    
                    IndexOut    = intersect(t(TriNum, 1:3)', index_out);    %    index into node(s) of TriNum that are outside                                  
                    if length(IndexOut)==2                  %    two nodes of TriNum are outside and one (target node) is inside    
                        OutNode1    = P(IndexOut(1), :);
                        OutNode2    = P(IndexOut(2), :);
                        InNode      = P(point, :);    
                        InterNode1 = (DIST(IndexOut(1))*InNode + DIST(point)*OutNode1)/(DIST(IndexOut(1)) + DIST(point));
                        InterNode2 = (DIST(IndexOut(2))*InNode + DIST(point)*OutNode2)/(DIST(IndexOut(2)) + DIST(point));
                        Numbers = Numbers + 2;
                        %   Construct subdivision triangles                        
                        tadd1 = [IndexOut(1)  Numbers-1  IndexOut(2)];
                        tadd2 = [IndexOut(2)  Numbers-1  Numbers-0];                    
                        tadd3 = [point        Numbers-1  Numbers-0];
                        NewTriangles    = [NewTriangles; tadd1; tadd2; tadd3];
                        NewNormals      = [NewNormals; normals(TriNum, :); normals(TriNum, :); normals(TriNum, :)];
                        Remove          = [Remove; TriNum];
                        NewNodes        = [NewNodes; InterNode1; InterNode2];  
                        if norm(InterNode1-InterNode2)<1024*eps
                            'here'
                        end
                    end                    
                    if (length(IndexIn)==1)&&(length(IndexOut)==1)          %   two nodes of TriNum are inside and one outside  
                        if isempty(intersect(Remove, TriNum))
                            InNode1 = P(IndexIn, :);
                            InNode2 = P(point, :);
                            OutNode = P(IndexOut, :);
                            InterNode1 = (DIST(IndexIn)*OutNode + DIST(IndexOut)*InNode1)/(DIST(IndexIn) + DIST(IndexOut));
                            InterNode2 = (DIST(point)*OutNode   + DIST(IndexOut)*InNode2)/(DIST(point) + DIST(IndexOut));
                            Numbers = Numbers + 2;
                            %   Construct subdivision triangles                
                            tadd1 = [IndexIn   Numbers-1  point];
                            tadd2 = [point     Numbers-1  Numbers-0];
                            tadd3 = [IndexOut  Numbers-1  Numbers-0]; % good
                            NewTriangles    = [NewTriangles; tadd1; tadd2; tadd3];   
                            NewNormals      = [NewNormals; normals(TriNum, :); normals(TriNum, :); normals(TriNum, :)];
                            Remove          = [Remove; TriNum];
                            NewNodes        = [NewNodes; InterNode1; InterNode2];                              
                        end
                    end                    
                end                
            end
        end
    end     
   
    %   Construct the new mesh
    t(Remove, :)        = [];
    normals(Remove, :)  = [];
    t                   = [t; NewTriangles];
    normals             = [normals; NewNormals];
    P                   = [P; NewNodes];
    
    %   Remove triangles with coincident points from the mesh (when two nodes are at the boundary)
    A = meshareas(P, t);
    index = find(A<1e-9);
    t(index, :)                 = [];
    normals(index, :)           = []; 
        
    %  Remove duplicated nodes from the mesh        
    [P, t]  = fixmesh(P, t);
    
    %   Remove duplicated triangles with equal centers
    C = meshtricenter(P, t);
    [~, index, ~] = unique(C, 'rows');
    t             = t(index, :);
    normals       = normals(index, :);
    
    %   Reorient triangles as required
    t = meshreorient(P, t, normals);
    
    %   Select electrodes  
    IndicatorElectrodes = zeros(size(t, 1), 1);    
    C = meshtricenter(P, t);
    for m =  1:NumberOfElectrodes
        %   Identify new electrode triangles     
        [DIST, NEARESTPOINT]    = meshimprint_distance(C, strge, m);
        IndicatorElectrodes(DIST<-0.1*avgedgelength) = m;    
    end   

    %  Finally, put electrodes up front sequentially (1, 2, 3, etc.)
    tt = [];
    nn = [];
    ie = [];
    for m = 1:strge.NumberOfElectrodes 
        index  = IndicatorElectrodes==m;
        tt     = [tt; t(index, :)];
        nn     = [nn; normals(index, :)];
        ie     = [ie; m*ones(sum(index), 1)];
    end
    index               = IndicatorElectrodes==0;
    t                   = [tt; t(index, :)];
    normals             = [nn; normals(index, :)];
    IndicatorElectrodes = [ie; IndicatorElectrodes(index)];
    
    index               = IndicatorElectrodes>0;
    t                   = [t(index, :); t(~index, :)];
    normals             = [normals(index, :); normals(~index, :)];
    IndicatorElectrodes = IndicatorElectrodes(index);
    IndicatorElectrodes(end+1:size(t, 1)) = 0;
end

function [si, vi] = vertices(P, t)
    %   si - triangles attached to every vertex (neighbor triangles)
    %   vi - Vertices attached to every vertex (neighbor edges)   
    si = cell(size(P, 1), 1);
    vi = cell(size(P, 1), 1);
    for m = 1:size(P, 1)
        temp  = (t(:, 1)==m)|(t(:, 2)==m)|(t(:, 3)==m);
        si{m} = find(temp>0);
        temp  = unique([t(si{m}, 1); t(si{m}, 2); t(si{m}, 3)]);
        temp(temp==m) = [];  
        vi{m} = temp(temp>0);
    end    
end

function [DIST, NEARESTPOINT] = meshimprint_distance(P, strge, Number)
%   Signed distance (DIST) from boundary for any point of the mesh
%   This distance is negative inside and positive outside
%   NEARESTPOINT is the nearest point at the boundary
%   Boundary nodes
    M = 100;
    R = strge.R;
    if strge.type(Number) == 1  %   Ring
        phi       = linspace(0, 2*pi, M);
        PB1(:, 1) = R*cos(phi);
        PB1(:, 2) = R*sin(phi);
        PB1(:, 3) = ones(M, 1)*strge.z1(Number);
        PB2(:, 1) = R*cos(phi);
        PB2(:, 2) = R*sin(phi);
        PB2(:, 3) = ones(M, 1)*strge.z2(Number);
        PB = [PB1; PB2];
        %   Find unsigned distance
        DT              = pdist2(PB, P);
        [DIST, point]   = min(DT);
        DIST            = DIST';
        %   Find signed distance
        index  = P(:, 3)>strge.z1(Number)&P(:, 3)<strge.z2(Number);
        DIST(index) = - DIST(index);
        %   Find nearest point at the boundary
        NEARESTPOINT = PB(point, :);
    end
    if strge.type(Number) == 2  %   Sector
        PB1(:, 1) = R*ones(M, 1)*cos(strge.phi1(Number));
        PB1(:, 2) = R*ones(M, 1)*sin(strge.phi1(Number));
        PB1(:, 3) = linspace(strge.z1(Number), strge.z2(Number), M)';
        PB2(:, 1) = R*ones(M, 1)*cos(strge.phi2(Number));
        PB2(:, 2) = R*ones(M, 1)*sin(strge.phi2(Number));
        PB2(:, 3) = linspace(strge.z1(Number), strge.z2(Number), M)';
        phi       = linspace(strge.phi1(Number), strge.phi2(Number), M);
        PB3(:, 1) = R*cos(phi);
        PB3(:, 2) = R*sin(phi);
        PB3(:, 3) = ones(M, 1)*strge.z1(Number);
        PB4(:, 1) = R*cos(phi);
        PB4(:, 2) = R*sin(phi);
        PB4(:, 3) = ones(M, 1)*strge.z2(Number);
        PB = [PB1; PB2; PB3; PB4];
        %   Find unsigned distance
        DT              = pdist2(PB, P);
        [DIST, point]   = min(DT);
        DIST            = DIST';
        %   Find signed distance
        angle       = atan(P(:, 2)./(P(:, 1)+eps));
        index       = P(:, 2)>0&P(:, 1)<0;
        angle(index)= angle(index) + pi;
        index       = P(:, 2)<0&P(:, 1)<0;
        angle(index)= angle(index) - pi;   
        index1  = angle>strge.phi1(Number)&angle<strge.phi2(Number);
        index2  = P(:, 3)>strge.z1(Number)&P(:, 3)<strge.z2(Number);
        index   = index1&index2;
        DIST(index) = - DIST(index);
        %   Find nearest point at the boundary
        NEARESTPOINT = PB(point, :);
    end
end