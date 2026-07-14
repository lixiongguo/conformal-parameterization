%% Declarations
phi=@(x) x.^3;
m=length(pointsOriginalData);
stretch=1.1;
kMatrix=zeros(2*m,2*n);


%% Create Random subset
    c=sort(randperm(length(pointsOriginalData),n));
    randomPoints=pointsOriginalData(c,:);  % output matrix
    randomNormals=normalsFromUser(c,:);
    
%% calculating the kMatrix 
for ii=1:m
    for jj=1:n
        kMatrix(ii,jj) = phi(norm(pointsOriginalData(ii,:)-randomPoints(jj,:)));
        kMatrix(ii,jj+n) = phi(norm(pointsOriginalData(ii,:)-(randomPoints(jj,:)+epsilon*randomNormals(jj,:))));
        kMatrix(ii+m,jj) = phi(norm((pointsOriginalData(ii,:)+epsilon*normalsFromUser(ii,:))-randomPoints(jj,:)));
        kMatrix(ii+m,jj+n) = phi(norm((pointsOriginalData(ii,:)+epsilon*normalsFromUser(ii,:))-(randomPoints(jj,:)+epsilon*randomNormals(jj,:))));
    end
end

%% calculating the weights vector
weights=kMatrix\[zeros(m,1); repmat(epsilon,m,1)];


%% Detect the size of the wanted containing cube and create the meshgrid.
distances=[abs(min(pointsOriginalData(:,1))) abs(max(pointsOriginalData(:,1))); 
           abs(min(pointsOriginalData(:,2))) abs(max(pointsOriginalData(:,2))); 
           abs(min(pointsOriginalData(:,3))) abs(max(pointsOriginalData(:,3)))];
cubeX=linspace(-stretch*distances(1,1),stretch*distances(1,2),gridSize);
cubeY=linspace(-stretch*distances(2,1),stretch*distances(2,2),gridSize);
cubeZ=linspace(-stretch*distances(3,1),stretch*distances(3,2),gridSize);
[x,y,z]=meshgrid(cubeX,cubeY,cubeZ);


%% create the [(p(i)+epsilon*n(i))] vector
p_i=[randomPoints; randomPoints+epsilon*randomNormals];

%% Calculate the actual Signed Distance Function for points in the meshgrid
SDF=zeros(gridSize,gridSize,gridSize);
for ii=1:2*n
    SDF=SDF + weights(ii)*(((x-p_i(ii,1)).^2+(y-p_i(ii,2)).^2+(z-p_i(ii,3)).^2).^(3/2));
end

[f,v]=isosurface(x,y,z,SDF,0);
numberOfFaces=length(f);
numberOfVertices=length(v);
fCount=size(f,2);


