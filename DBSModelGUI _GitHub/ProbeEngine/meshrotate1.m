function PP = meshrotate1(P, nx, ny, nz)
%   Rotation (Rodrigues' rotation formula)
%   If the original structure has the normal vector (or the axis) in the
%   z-direction; it is so rotated that the normal vector now becomes [nx ny nz]

%   Copyright SNM 2017-2020

    temp    = sqrt(nx^2 + ny^2 + nz^2);
    nx      = nx/temp;
    ny      = ny/temp;
    nz      = nz/temp;
    
    K0    = repmat([0 0 0]', [1 size(P, 1)])';  
    theta = acos(nz);
    k     = [-ny +nx 0]';
    K     = K0;
    if dot(k, k)>1e-6
        K     = repmat(k, [1 size(P, 1)])'/sqrt(dot(k, k));
    end
    PP = P*cos(theta) + cross(K, P, 2)*sin(theta) + K.*repmat(dot(K, P, 2), [1 3])*(1-cos(theta));
end