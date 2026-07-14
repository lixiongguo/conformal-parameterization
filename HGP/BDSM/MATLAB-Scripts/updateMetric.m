function [newMetric] = updateMetric(mesh_EL, mesh_EV, phi) 
% % This function takes each edge, finds it's connected vertices, and calculates
% % the value of newLength = oldLength * exp((phi(v1) + phi(v2))/2).
% % By discretization, this is the new metric.

phi_i = exp(phi(mesh_EV(:,1))/2);
phi_j = exp(phi(mesh_EV(:,2))/2);

newMetric=((phi_i.*phi_j)'*sparse(1:length(mesh_EL), 1:length(mesh_EL), mesh_EL))';

