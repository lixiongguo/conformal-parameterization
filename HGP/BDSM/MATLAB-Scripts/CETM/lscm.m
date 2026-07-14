function x = lscm(el, t)

nv = max(max(t));
nf = size(t, 1);

% MA = sparse( t(:, [1 2 3 2 3 1]), nv+t(:, [2 3 1 1 2 3]), [ones(nf, 3) -ones(nf, 3)], nv*2, nv*2 );
MA = sparse( [t t+nv], [nv+t(:, [2 3 1]) t(:,[2 3 1])], [ones(nf, 3) -ones(nf, 3)], nv*2, nv*2 );
MA = (MA+MA')/2;

L = cotLaplaceFromFaceAngles( meshAnglesFromFaceEdgeLen2(el.^2), t, nv );

%% conformal laplacian
CL = blkdiag( L, L ) - MA*2;

CL(t(1,[1 1 2 2])+[0 nv 0 nv], :) = sparse( 1:4, t(1,[1 1 2 2])+[0 nv 0 nv], ones(1,4), 4, nv*2 );
b = zeros(nv*2, 1);
b( t(1,2)+nv ) = el(1,3);

x = CL\b;
x = reshape( x, [], 2 );
