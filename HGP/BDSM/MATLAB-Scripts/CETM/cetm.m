function [UV, output_edge_lengths] = cetm(mesh_V,mesh_F, mesh_FE, ELS, maxNumIter, kappaIndices, kappa)

nv = size(mesh_V, 1);
nf = size(mesh_F, 1);


VV2T = sparse(mesh_F, mesh_F(:, [2 3 1]), repmat(1:nf, 1, 3), nv, nv);
VV2T( VV2T & VV2T' ) = 0;
[B, ~] = find(VV2T);


EL = mesh_FE;
l = arrayfun(@(y) ELS(y), EL);

Theta = ones(nv, 1)*2*pi; %set prescribed angles for all vertices (internal and boundary) to 2pi
ufixed = [B zeros(numel(B),1)]; %fix the conformal scale factor of the boundary vertices to zero
notB = [1:nv];
notB(B) = [];
if(nargin >= 7)
     Theta(B) = pi - target_curvature(B);
     Theta(nB) = Theta(nB) - target_curvature(nB);
    %obtain the indices of boundary vertices for which we don't specify curvature
    scaleFixIndices = B;
    ufixed = [scaleFixIndices zeros(numel(scaleFixIndices),1)];
end



[UV, nbrokentris, lcrerr, flaterr, output_edge_lengths] = dcflatten_wrap(mesh_FV, nv, l, 100, 1, Theta,ufixed, 0);
% figure;
% trimesh(mesh_F, mesh_V(:,1), mesh_V(:,2));
% figure;
% trimesh(mesh_F, UV(:,1), UV(:,2));

end