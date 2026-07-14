function [out] = findDiffNew(a,b,c,d,e,f,vertex, mesh_EV)
%FINDDIFF Summary of this function goes here
%   Detailed explanation goes here
out = [a b c d e f];
g = mode(out);
% Remain with only the uncommon edges
out(out == g) = [];

% Rearrange the first 2 remaining edges so that the first constains the
% first vertex of the common edge.
if (ismember(vertex, mesh_EV(out(2),:)));
    [out(1:2)] = swap(out(1),out(2));
end
if (ismember(vertex, mesh_EV(out(4),:)));
    [out(3:4)] = swap(out(3),out(4));
end

out = [out g];

