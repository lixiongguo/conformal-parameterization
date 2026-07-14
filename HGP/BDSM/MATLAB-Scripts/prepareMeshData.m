function [ mesh_ConnectedTriangles, mesh_EF, mesh_ELS, ...
    mesh_EV, mesh_FacesAdjacency, mesh_FacesAreas, mesh_FE, mesh_FV, ...
    mesh_InnerEdges, mesh_V, mesh_VertexArea, mesh_VerticesAdjacency, ...
    meshName, Test_OrientationResults] = prepareMeshData()
%%
% set up field types
v = []; vt = []; vn = []; f.v = []; f.vt = []; f.vn = [];
cones = [];

fid = fopen(fname);

% parse .obj file
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end  % exit at end of file
    ln = sscanf(tline,'%s',1); % line type
    %disp(ln)
    switch ln
        case 'v'   % mesh vertexs
            v = [v; sscanf(tline(2:end),'%f')'];
        case 'vt'  % texture coordinate
            vt = [vt; sscanf(tline(3:end),'%f')'];
        case 'vn'  % normal coordinate
            vn = [vn; sscanf(tline(3:end),'%f')'];
        case 'f'   % face definition
            fv = []; fvt = []; fvn = [];
            str = textscan(tline(2:end),'%s'); str = str{1};
            
            nf = length(findstr(str{1},'/')); % number of fields with this face vertices
            
            [tok str] = strtok(str,'//');     % vertex only
            for k = 1:length(tok) fv = [fv str2num(tok{k})]; end
            
            if (nf > 0)
                [tok str] = strtok(str,'//');   % add texture coordinates
                for k = 1:length(tok) fvt = [fvt str2num(tok{k})]; end
            end
            if (nf > 1)
                [tok str] = strtok(str,'//');   % add normal coordinates
                for k = 1:length(tok) fvn = [fvn str2num(tok{k})]; end
            end
            f.v = [f.v; fv]; f.vt = [f.vt; fvt]; f.vn = [f.vn; fvn];
        case 'c'
            c = sscanf(tline(2:end),'%f %f');
            c(2) = c(2) * 2; % now it's equal to the cone valence
            cones = [cones; c'];
    end
end
fclose(fid);

% set up matlab object
obj.v = v; obj.vt = vt; obj.vn = vn; obj.f = f;
obj.cones = cones;
nv = size(obj.v,1);
nf = size(obj.f.v,1)




% Prepare the variables expected in boundedDistortion.m
mesh_V1 = obj.v;
mesh_FV1 = f.v;
mesh_FE1 = zeros(size(f.v,1),3);
mesh_EV1 = zeros(size(f.v,1),2);

mesh_VerticesAdjacency1 = sparse(size(obj.v,1),size(obj.v,1));
count = 1;

for ii=1:size(f.v,1)
    v = f.v(ii,:);
    if (mesh_VerticesAdjacency1(v(1), v(2)) == 0)
        mesh_VerticesAdjacency1(v(1), v(2)) = count;
        mesh_VerticesAdjacency1(v(2), v(1)) = count;
        mesh_EV1(count,:) = [v(1) v(2)];
        count = count + 1;
    end
    if (mesh_VerticesAdjacency1(v(2), v(3)) == 0)
        mesh_VerticesAdjacency1(v(2), v(3)) = count;
        mesh_VerticesAdjacency1(v(3), v(2)) = count;
        mesh_EV1(count,:) = [v(2) v(3)];
        count = count + 1;
    end
    if (mesh_VerticesAdjacency1(v(3), v(1)) == 0)
        mesh_VerticesAdjacency1(v(3), v(1)) = count;
        mesh_VerticesAdjacency1(v(1), v(3)) = count;
        mesh_EV1(count,:) = [v(3) v(1)];
        count = count + 1;
    end
    mesh_FE1(ii,1) = mesh_VerticesAdjacency1(v(2), v(3));
    mesh_FE1(ii,2) = mesh_VerticesAdjacency1(v(1), v(3));
    mesh_FE1(ii,3) = mesh_VerticesAdjacency1(v(1), v(2));
end

ne = count;
mesh_InnerEdges = ones(ne,1);


VV2T = sparse(mesh_FV1, mesh_FV1(:, [2 3 1]), repmat(1:nf, 1, 3), nv, nv);
VV2T( VV2T & VV2T' ) = 0;
[B, ~] = find(VV2T);
mesh_VertexArea1 = 2 * pi * ones(size(obj.v,1),1);
mesh_VertexArea1(B) = pi;
a = nchoosek(B,2)

for ii=1:size(a,1)
    if(mesh_VerticesAdjacency1(a(ii,1),a(ii,2)))
        mesh_InnerEdges(mesh_VerticesAdjacency1(a(ii,1),a(ii,2))) = 0;
    end    
end


mesh_EF1 = sparse([mesh_FE1(:,1); mesh_FE1(:,2); mesh_FE1(:,3)], [1:nf ; 1:nf ; 1:nf],ones(nf*3,1),ne,nf);
mesh_EF1(mesh_EF1==0)=[];
%% POC
isequal(mesh_V1,mesh_V)
isequal(mesh_FV1,mesh_FV);
isequal(mesh_FE1,mesh_FE);
isequal(mesh_EV1,mesh_EV);
isequal(mesh_VertexArea1,mesh_VertexArea);
isequal(mesh_FE1,mesh_FE);

end 