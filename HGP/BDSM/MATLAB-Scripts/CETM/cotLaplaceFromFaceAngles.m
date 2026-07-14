function L = cotLaplaceFromFaceAngles(angles, t, nv)

%angles(angles == 0) = 0.1;

% L = - sparse( t(:,[2 3 1]), t(:,[3 1 2]), cot(angles), nv, nv );
% L = L+L';
L = - sparse( t(:,[2 3 1 3 1 2]), t(:,[3 1 2 2 3 1]), cot([angles angles]), nv, nv );
L = spdiags(-sum(L,2), 0, L);

%Ofir's addition - not sure if this is the best way to go but at least it avoids complete failure due to inf when triangle inequality doesn't hold
L(~isfinite(L)) = 0;