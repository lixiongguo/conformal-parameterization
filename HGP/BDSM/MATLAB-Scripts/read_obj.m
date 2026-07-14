function [target_curvature] = readObj(pathToOpen, meshName, pathToSave)
fname = 'C:\Users\Nadav\Desktop\HGP_Results\block.obj';
% obj = readObj(fname)


% set up field types
v = 0; cones = [];

fid = fopen(pathToOpen);

% parse .obj file
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end  % exit at end of file
    ln = sscanf(tline,'%s',1); % line type
    %disp(ln)
    switch ln
        case 'v'   % mesh vertexs
            v = v + 1;;
        case 'c'
            c = sscanf(tline(2:end),'%f %f');
            cones = [cones; c'];
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

target_curvature = sparse(v,1);
target_curvature(cones(:,1)) = 2*pi - cones(:,2) * pi;

filename = strcat(pathToSave, meshName, '.mat');
save(filename, 'target_curvature', 'meshName');
clear all; 



