function [ELS_next, candidate_norm] = setStartConditions(ELS_current, step_candidate, target_curvature, current_norm, nearZero, mesh_FE, mesh_FV, mesh_V, mesh_VertexArea, interiorVertices)
%SETSTARTCONDITINONS Summary of this function goes here
%   Detailed explanation goes here
    candidate_curvature = calculateGaussianCurvatureBoundedDistortion(ELS_current + step_candidate, mesh_FE, mesh_FV, mesh_V, mesh_VertexArea, nearZero);
    %candidate_norm = max(abs(candidate_curvature(interiorVertices) - target_curvature(interiorVertices)));
    candidate_norm = norm(candidate_curvature(interiorVertices)- target_curvature(interiorVertices),Inf);
    if (((abs(candidate_norm) - abs(current_norm))<nearZero) || isequal(step_candidate,zeros(length(step_candidate),1)))
        ELS_current(abs(ELS_current) < nearZero) = 0;
        step_candidate(abs(step_candidate) < nearZero) = 0;
        ELS_next = ELS_current + step_candidate;
    else
        step_candidate(abs(step_candidate) < nearZero) = 0;
        [ELS_next, candidate_norm]  = setStartConditions(ELS_current, step_candidate .* 0.5, target_curvature, current_norm, nearZero, mesh_FE, mesh_FV, mesh_V, mesh_VertexArea, interiorVertices);
    end
end

