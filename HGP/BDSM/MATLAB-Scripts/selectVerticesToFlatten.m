function [K_new] = selectVerticesToFlatten(calculatedGaussianCurvature, epsilon)
K_new = 0;
%Add index to the curvatues vector so it is easily restored later.
K = [(1:length(calculatedGaussianCurvature))' calculatedGaussianCurvature];

%Set numbers - constrained and unconstrained vertices
constVertices = floor(epsilon*length(K));
unconstVertices = length(K) - constVertices;

fprintf('const vertices is %i\n',constVertices);

%Make sure that the epsilon isn't actually 0.
if (constVertices<1)
    return
end

% sort the curvature vector so we can easily choose the minimal onse.
K_sorted = sortrows(abs(K), 2);

%Now use the original K to preserve the sign
K=K(K_sorted(:,1),:);

%Extract the value of the curvature to be added to each unconstrained
%vertex
curvatureToRedistribute = sum(K(1:constVertices,2))./unconstVertices;

%Preparing the redistribute matrix
vectorToRedistribute = [zeros(unconstVertices, 1) (curvatureToRedistribute * ones(unconstVertices, 1))];

%Redistribute curvature to the wanted vertices
K_new = [K(1:constVertices,:); K(constVertices+1:end,:) + vectorToRedistribute];

%zero the curvature for smaller than epsilon vertices
K_new(1:constVertices,2) = zeros;

%prepare the vector for output
K_new = sortrows(K_new, 1);
K_new = K_new(:,2);
return
