function [ commonEdges ] = nonTreeEdge2meshEdge(mesh_FacesAdjacency, mesh_FE)
% This function receives a list of edges that are not included in the
% spanning tree, and returns the indices of their corrolating edges in the
% original mesh

adjacencyGraph = mesh_FacesAdjacency;
adjacencyGraph(adjacencyGraph~=0) = 1;
ANS = graphminspantree(adjacencyGraph);
cutCandidate = adjacencyGraph-ANS-ANS';
cutCandidate = cutCandidate + cutCandidate';
cutCandidate(cutCandidate~=0) = 1;
commonEdges = unique(mesh_FacesAdjacency(find(cutCandidate)));

end

