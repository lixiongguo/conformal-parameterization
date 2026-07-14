%% --------------------------------------------------------------------%
%                                                                      %
%                            Initialization                            %
%                                                                      %
%----------------------------------------------------------------------%
% Pre determining the solver type.
% cvx_solver sedumi
% cvx_solver mosek

% Initializing variables
tic;
maxColorSets = 20;
K = zeros(length(mesh_V), 1);
phi = zeros(length(mesh_V), 1);
epsilon = pi/200;
threshold = 0;
nearZero = 10^(-7);
ii = 1;
jj = 1;
interior_vertices = find(mesh_VertexArea == 2*pi);
L_interior = mesh_L(interior_vertices, :);
coneCandidates = zeros(length(K(:,1)),1);
flatAndUnflatVertices = zeros(length(mesh_V),1);

% Calculate initial Gaussian curvature
calculated_VA = calculateGaussianCurvature(mesh_EL(:,1), mesh_FE, mesh_F, mesh_V);
K(:,1) = mesh_VertexArea - calculated_VA;

% get rid of values that may be numerical errors, set by parameter nearZero
K(abs(K(:,1))<nearZero) = 0;
numOfCones(1) = nnz(K(:,1));
flatAndUnflatVertices (K == 0) = 1; 

nextIter = true;
coneOrNotCone = zeros(length(mesh_V),1);

%% --------------------------------------------------------------------%
%                                                                      %
%                                Flattening                            %
%                                                                      %
%----------------------------------------------------------------------%


while ( nextIter == true && (threshold <=(pi/4)))

    % Increase the threshold
    threshold = threshold + epsilon;
        
    % Prepare a list of unflat vertices to see if seperate
    [nextIter, coneOrNotCone] = conesAreSeperated(flatAndUnflatVertices(:,ii), coneOrNotCone, mesh_F, ii, interior_vertices);
    if (nextIter == 0)
        flatVertices = find(~K(:,ii));
        numOfCones(ii) = nnz(K(:,ii));
        break;
    end
    
    % Create a vector of vertices to be falttened.
    flatVertices = find(abs(K(:,ii)) < threshold);
    flatVertices = flatVertices(~ismember(flatVertices,find(coneOrNotCone)));
    
    %Calculate the optimization problem based on the constraints.
    cvx_begin quiet
        variable cvx_phi(length(mesh_V))
        variable k_new(length(mesh_V))
        minimize norm(mesh_Areas*cvx_phi, 2)
        subject to
            L_interior*cvx_phi - k_new(interior_vertices) == - K(interior_vertices,1);
            k_new(flatVertices) == 0;
    cvx_end    

    phi(:,ii) = cvx_phi;

    % fix numerical inaccuracies
    k_new(abs(k_new)<nearZero) = 0;
    K(:,ii+1) = k_new;
    
    %Prepare next iteration
    flatAndUnflatVertices (K(:,ii+1)==0,ii+1) = 1;
    ii = ii + 1;
    numOfCones(ii) = nnz(K(:,ii));
 end

   
    
%% --------------------------------------------------------------------%
%                                                                      %
%                                Rounding    	                       %
%                                                                      %
%----------------------------------------------------------------------%
    

    
    conesToRound = find(K(:,ii));
    conesToRound(:,2) = mod(K(conesToRound(:,1),ii),(pi/2)); %sort by closest to k*(pi/2)
    conesToRound = sortrows(conesToRound,2);
    roundedCones=zeros(length(conesToRound),1);
    for jj = 1:size(conesToRound,1)
         roundedCones(jj) = round(K(conesToRound(jj,1) ,ii+jj-1)/(pi/2))*pi/2; %ii+jj-1
         if  (K(conesToRound(jj,1) ,ii+jj-1) / (pi/2) == round(K(conesToRound(jj,1) ,ii+jj-1) / (pi/2)))
            
            phi(:,end + 1) = phi(:,end);
            K(:,end + 1) = K(:,end);
            flatAndUnflatVertices(:,end + 1) = flatAndUnflatVertices(:,end);
            numOfCones(end + 1) = numOfCones(end);
             continue;
         end
         
         cvx_begin quiet
            variable cvx_phi(length(mesh_V))
            variable k_new(length(mesh_V))
            minimize norm(mesh_Areas*cvx_phi, 2)
            subject to
                L_interior*cvx_phi - k_new(interior_vertices) == - K(interior_vertices,1);
                k_new(flatVertices) == 0;
                k_new(conesToRound(1:jj,1)) == roundedCones(1:jj);
        cvx_end  
        phi(:,ii+jj) = cvx_phi;
        
        % fix numerical inaccuracies
        k_new(abs(k_new)<nearZero) = 0;
        K(:,ii+jj) = k_new;
        flatAndUnflatVertices (K(:,ii+jj)==0,ii+jj) = 1;
        numOfCones(ii+jj) = nnz(K(:,ii+jj));
    end
target_curvature = K(:,end);

%clearvars -except target_curvature




%% --------------------------------------------------------------------%
%                                                                      %
%                   prepare data to be sent back to Maya               %
%                                                                      %
%----------------------------------------------------------------------%

if (ii > maxColorSets)
    export = unique([1 round(size(K(:,1:ii+jj),2)/(maxColorSets)):round(size(K(:,1:ii+jj),2)/(maxColorSets)):size(K(:,1:ii+jj),2) size(K,2)]);
    numColorSets = length(export);
    K_to_MAYA = K(:, export);
    numOfCones_to_MAYA = numOfCones(export);
    flatAndUnflatVertices_to_MAYA = flatAndUnflatVertices(:, export);
    maxCurvature=max(K(:, export));
    minCurvature=min(K(:, export));
else
    export = 1:size(K,2);
    numColorSets = ii+1;
    K_to_MAYA = K;
    numOfCones_to_MAYA = numOfCones;
    flatAndUnflatVertices_to_MAYA = flatAndUnflatVertices;
    maxCurvature = max(K);
    minCurvature = min(K);
end
numOfIterations = ii+jj;

%% --------------------------------------------------------------------%
%                                                                      %
%                          Proof of Correctness                        %
%                                                                      %
%----------------------------------------------------------------------%
figure (1)
graphs = gcf;
set(gcf, 'Position', get(0, 'Screensize'));

subplot(2,2,3);       
newMetric = updateMetric(mesh_EL, mesh_EV, phi(:,end));
newCurvature = mesh_VertexArea - calculateGaussianCurvature(newMetric, mesh_FE, mesh_F, mesh_V);
plot(-newCurvature(interior_vertices));
grid on;
title('Computed new curvature, based on evolved phi to update metric');
xlabel('Index of vertex');
ylabel('Curvature value');
ax=gca;


subplot(2,2,1);

plot (K(interior_vertices,1));
xlabel('Index of vertex');
ylabel('Curvature value');
hold on;
grid on;
plot (K(interior_vertices,end));
title('Incremental Flattening curvature values, final result');
xlabel('Index of vertex');
ylabel('Curvature value');
grid on;
legend('Initial Iteration', 'Final Iteration');
hold off;


subplot(2,2,2);   
IncFlat_dist=phi(:,end).^2;
IncFlat_dist(IncFlat_dist<nearZero) = 0;
hist(IncFlat_dist);
title('Incremental Flattening Distortion histogram, phi^2');
xlabel('Value of phi^2');
xlabel('Number of vertices');


subplot(2,2,4);       
ARAP_dist = (exp(phi(:,end))-1).^2;
hist(ARAP_dist);
title('ARAP Distortion histogram, (exp(phi)-1).^2');
xlabel('Value of phi^2');
xlabel('Number of vertices');
runtime=toc;

%% Section for batch processing InfFlat

pathToSaveFig = strcat('C:\code\SpaceDeformation2D\results\Incremental-Flattening-Figures\', meshName, '.jpg')
saveas(graphs, pathToSaveFig);
clearvars -except target_curvature meshName
pathToSaveVars = strcat('C:\code\SpaceDeformation2D\results\target_curvature\', meshName,'-target-curvature.mat');
save(pathToSaveVars);
close all;

