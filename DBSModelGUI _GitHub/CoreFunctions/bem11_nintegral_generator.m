function bem11out = bem11_nintegral_generator(bem02out)

    % From structure to individual fields
    cellfun(@(f) assignin('caller', f, bem02out.(f)), fieldnames(bem02out));

    %%   Add accurate integration for electric field/electric potential on neighbor facets
    %   Indexes into neighbor triangles
    RnumberE        = 32;      %   number of neighbor triangles for analytical integration of electric field
    RnumberP        = 32;      %   number of neighbor triangles for analytical integration of electric potential
    ineighborE      = knnsearch(Center, Center, 'k', RnumberE);   % [1:N, 1:RnumberE]
    ineighborP      = knnsearch(Center, Center, 'k', RnumberP);   % [1:N, 1:RnumberP]
    ineighborE      = ineighborE';          %   do transpose  
    ineighborP      = ineighborP';          %   do transpose  
    
    %%   Add accurate integration for electric field/electric potential on neighbor facets
    %   Indexes into neighbor triangles
    [PC, integralpd]    = meshneighborints_Pn(P, t, normals, Area, Center, RnumberP, ineighborP);
    EC                  = meshneighborints_En(P, t, normals, Area, Center, RnumberE, ineighborE);
    ec                  = isnan(EC); 
    indicator           = sum(sum(ec))
    EC(ec)              = 0;
    pc                  = isnan(PC); 
    indicator           = sum(sum(pc))
    PC(pc)              = 0;
    
    %%   Normalize sparse matrix EC by variable contrast (for speed up)
    N   = size(Center, 1);
    ii  = ineighborE;
    jj  = repmat(1:N, RnumberE, 1); 
    temp = Contrast(ineighborE);
    CO  = sparse(ii, jj, temp);
    EC  = CO.*EC;
    
    %%  Electrode preconditioner M (left)
    tic
    Ne          = sum(IndicatorElectrodes>0); 
    tempC       = Center(1:Ne, :); 
    tempA       = Area(1:Ne); 
    A           = repmat(tempA, 1, length(tempA));
    M           = (1/(4*pi))*1./dist(tempC').*A';       %   base preconditioner matrix
    for m = 1:Ne                                        %   base matrix with zero elements
        M(m, m) = 0;
    end
    for m = 1:Ne                                        %   put in neighbor integrals
        mneighbors          = ineighborE(:, m);
        mneighborsn         = mneighbors(mneighbors<=Ne);
        M(mneighborsn, m)   = M(mneighborsn, m) + integralpd(mneighbors<=Ne, m)/(4*pi);
    end
    M = inv(M);                                        %   direct inversion - replace
    disp([newline 'Preconditioner computed in ' num2str(toc) ' s']);

    % From individual fields to structure
    bem11out = struct('EC', EC, 'PC', PC, 'integralpd', integralpd, 'M', M);
end
