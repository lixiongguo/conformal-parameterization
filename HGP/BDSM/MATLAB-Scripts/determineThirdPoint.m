function [ UV UVindex ] = determineThirdPoint(mesh_EV, mesh_FE, mesh_FV, UV, UVindex, path, distSorted, ELS, mesh_VerticesAdjacency, nearZero, root)
    path{root}=[root root]; %so it takes (end-1)...
    sortedPath = [cellfun(@(y) y(end-1), path)' , (1:length(mesh_FE))'];
    sortedPath = sortedPath(distSorted(2,:),:);

    C2 = zeros(length(sortedPath),1);
    D2 = zeros(length(sortedPath),1);

    for ii=1:length(sortedPath)
       currentFace = sortedPath(ii, 2);
       lastFace = sortedPath(ii, 1);

        % find the third vertex of the current triangle and its possition in
        % the triplet
        %%

        commonVertices = mesh_EV(mode([mesh_FE(lastFace,:) mesh_FE(currentFace,:)]),:);
        if (ii==1)
            commonVertices = [mesh_FV(currentFace,1) mesh_FV(currentFace,2)];
        end

        thirdVertexOfCurrent  = setdiff(mesh_FV(currentFace,:),commonVertices);
        thirdVertexOfLast = setdiff(mesh_FV(lastFace,:),commonVertices);

        %find the place of each common vertex within UV triplet
        [ ~ , locationInPrev1] = ismember(commonVertices(1) , UVindex(((lastFace-1) * 3 + 1):((lastFace-1) * 3 + 3)));
        [ ~ , locationInPrev2] = ismember(commonVertices(2) , UVindex(((lastFace-1) * 3 + 1):((lastFace-1) * 3 + 3)));
        locationInPrev3 = 6 - locationInPrev1 - locationInPrev2;

        [ ~ , locationInOrigin1] = ismember(commonVertices(1) ,mesh_FV(currentFace,:));
        [ ~ , locationInOrigin2] = ismember(commonVertices(2) ,mesh_FV(currentFace,:));
        locationInOrigin3 = 6 - locationInOrigin1 - locationInOrigin2;

        %Copy the common veritices' uv values to the new face-UV triplet
        UV((currentFace-1) * 3 + locationInOrigin1, :) =  UV((lastFace - 1) * 3 + locationInPrev1, :);
        UV((currentFace-1) * 3 + locationInOrigin2, :) =  UV((lastFace - 1) * 3 + locationInPrev2, :);
        UVindex((currentFace - 1) * 3 + locationInOrigin1) = commonVertices(1);
        UVindex((currentFace - 1) * 3 + locationInOrigin2) = commonVertices(2);

        % Calculate the position of the uncommon vertex

        x0 = UV((currentFace-1) * 3 + locationInOrigin1, 1);
        y0 = UV((currentFace-1) * 3 + locationInOrigin1, 2);

        x1 = UV((currentFace-1) * 3 + locationInOrigin2, 1);
        y1 = UV((currentFace-1) * 3 + locationInOrigin2, 2);

        x_prev = UV((lastFace-1) * 3 + locationInPrev3, 1);
        y_prev = UV((lastFace-1) * 3 + locationInPrev3, 2);

        %%

        a_squared = ELS(mesh_VerticesAdjacency(commonVertices(1),thirdVertexOfCurrent));
        b_squared = ELS(mesh_VerticesAdjacency(commonVertices(2),thirdVertexOfCurrent));
        c_squared = ELS(mesh_VerticesAdjacency(commonVertices(1),commonVertices(2)));

        r = b_squared + c_squared - a_squared;

        C2(ii) = ((r/((2 * c_squared)))+ i*sqrt((b_squared/c_squared)-((r^2)/(4 * (c_squared^2))))) * ((x0 - x1) + i*(y0 - y1)) + (x1 + i*y1);
        D2(ii) = ((r/((2 * c_squared)))- i*sqrt((b_squared/c_squared)-((r^2)/(4 * (c_squared^2))))) * ((x0 - x1) + i*(y0 - y1)) + (x1 + i*y1);
        dist = sqrt(((x_prev - real(C2(ii)))^2)+((y_prev - imag(C2(ii)))^2));
        dist_op = sqrt(((x_prev - real(D2(ii)))^2)+((y_prev - imag(D2(ii)))^2));

        if(dist >= dist_op)
           x = real(C2(ii));
           y = imag(C2(ii));
        else
           x = real(D2(ii));
           y = imag(D2(ii));
        end
         % assign the new values to UV and UVindex
        UV((currentFace-1) * 3 + locationInOrigin3, 1) =  x;
        UV((currentFace-1) * 3 + locationInOrigin3, 2) =  y;
        UVindex((currentFace-1) * 3 + locationInOrigin3) =  thirdVertexOfCurrent;

        % Maya will take a sequetial UVindex vector....

    end
    UVindex = 1:size(UV,1);
end

