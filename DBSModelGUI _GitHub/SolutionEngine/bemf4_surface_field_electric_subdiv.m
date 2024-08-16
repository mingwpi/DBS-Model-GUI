function [P, E] = bemf4_surface_field_electric_subdiv(c, P, t, Area, mode, modeArg) 
    %% Input parsing
    if(nargin < 4)
        error('Not enough input arguments');
    end
    %Assign default mode
    if(nargin < 5)
        mode = 'gauss';
    end
    
    %Assign default subdivision parameter
    if(nargin < 6)
        if(strcmp(mode, 'gauss'))
            modeArg = 7;
        else
            modeArg = 3;
        end
    end
    
    %Assign actual internal parameters
    if(strcmp(mode, 'gauss'))
        gauss = modeArg;
    else
        gauss = 0;
        subdivParam = modeArg;
    end
    
    %% Compute model subdivision parameters
                        %   number of integration points in the Gaussian quadrature  
                        %   for the outer potential integrals
                        %   Numbers 1, 4, 7, 13, 25 are permitted 
    %   Gaussian weights for analytical integration (for the outer integral)
    if      gauss == 0  
        [coeffS, weightsS, IndexS]  = tri(subdivParam);
    elseif  gauss == 1  
        [coeffS, weightsS, IndexS]  = tri(1, 1);
    elseif  gauss == 4  
        [coeffS, weightsS, IndexS]  = tri(4, 3);
    elseif  gauss == 7
        [coeffS, weightsS, IndexS]  = tri(7, 5);
    elseif  gauss == 13 
        [coeffS, weightsS, IndexS]  = tri(13, 7);
    elseif  gauss == 25 
        [coeffS, weightsS, IndexS]  = tri(25, 10);
    else
        error('Invalid Gaussian subdivision parameter');
    end
    
    % Subdivide all model triangles according to subdivision parameters
    Center_subdiv = zeros(IndexS*size(t, 1), 3);
    c_subdiv = zeros(IndexS*size(c, 1), 1);
    Area_subdiv = zeros(IndexS*size(Area, 1), 1);
    P1 = P(t(:,1), :);
    P2 = P(t(:,2), :);
    P3 = P(t(:,3), :);
    for j = 1:IndexS
        currentIndices = ([1:size(t, 1)] - 1)*IndexS + j; 
        Center_subdiv(currentIndices, :) = coeffS(1, j)*P1 + coeffS(2, j)*P2 + coeffS(3, j)*P3;
        c_subdiv(currentIndices, :) = c;    %Charge density doesn't change!
        Area_subdiv(currentIndices, :) = weightsS(j)*Area;
    end

    [P, E] = bemf4_surface_field_electric_plain(c_subdiv, Center_subdiv, Area_subdiv);

    %Every column contains the subdivided quantities for one full triangle
    P_temp = reshape(P, IndexS, []);
    Ex_temp = reshape(E(:,1), IndexS, []);
    Ey_temp = reshape(E(:,2), IndexS, []);
    Ez_temp = reshape(E(:,3), IndexS, []);

    %Eavg = integral(E dA)/A.  dA = subdivided area. subdivided area/A = weightsS
    P = transpose(weightsS*P_temp);
    Ex_avg = transpose(weightsS*Ex_temp);
    Ey_avg = transpose(weightsS*Ey_temp);
    Ez_avg = transpose(weightsS*Ez_temp);

    E = [Ex_avg, Ey_avg, Ez_avg];
end