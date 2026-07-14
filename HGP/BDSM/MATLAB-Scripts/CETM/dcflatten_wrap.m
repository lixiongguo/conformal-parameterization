function [y, nbrokentris, lcrerr, flaterr, output_edge_lengths, output]= dcflatten_wrap(t, nv, edgelens, maxNumIter, fid, varargin)

% function [y, nbrokentris, lcrerr]= dcflatten_wrap(t, nv, edgelens, Theta, ufix)
%
% wrapper for dcflatten
%
% Input parameters:
%
% t - mesh connectivity (list of triplets)
% nv - number of vertices
% edgelens - length of edges in each triangle, given mesh (x,t), can be computed by sqrt( meshFaceEdgeLen2s(x, t) )
% maxNumIter - maximum number of iterations
% fid - file id for output. use 1 for standart output or use fopen/fclose to write to a file
% Theta - angle sum for each vertex
% ufix - indices and value of u to be fixed
% bLSCM - layout with LSCM
%
% for natural boundary condition, Theta and ufix can be ignored

%if exist('mosekopt')
%    warning('Mosek in path, will interfere with fminunc, remove it from path first!');
%end

[y, nbrokentris, flaterr, output_edge_lengths, output] = dcflatten(t, nv, edgelens, maxNumIter, fid, varargin{:});

if nargout>2
    %% sanity check
    llcr_input = LCRFromFaceEdgeLens(t, nv, edgelens);
    llcr_output = LCRFromFaceEdgeLens(t, nv, output_edge_lengths);
    llcr_output_embedding = LCRFromFaceEdgeLens(t, nv, sqrt(meshFaceEdgeLen2s(y, t)));

    if 0
        % won't work on closed meshes
        lcrerr = norm(llcr_input - llcr_output_embedding, 'inf'); %compare the original cross ratios to the cross ratios of the final embedding
    else
        lcrerr = norm(llcr_input - llcr_output, 'inf'); %compare the original cross ratios to the cross ratios of lengths obtained from the convex optimization (before embedding)
    end
    if lcrerr > 1e-3
        warning( ['length cross ratios are not reproduced by CETM, maximum error is: ' num2str(lcrerr)] );
    end
end

end


 
function llcr = LCRFromFaceEdgeLens(t, nv, edgelens)
%length cross ratios for each edge.
%the output is given by an nv X nv sparse adjacency matrix.
%boundary edges (both halfedges) are set to 0 since there is no cross ratio defined for them.

halfedges = sparse(t, t(:,[2 3 1]), 1, nv, nv);
boundaryEdges = xor(halfedges, halfedges');
llcr = sparse(t, t(:,[2 3 1]), reallog(edgelens(:,[2 3 1])./edgelens), nv, nv); %for each halfedge compute the corresponding log length ratio
llcr = llcr + llcr'; %sum the two parts of the log length ratios to obtain the cross ratios
llcr(boundaryEdges) = 0; %set the cross ratio on boundary edges (both halfedges) to 0 since it is not well defined there
end
