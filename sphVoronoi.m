function [voronoi, duplicates] = sphVoronoi(dirs, faces)
%SPVORONOI Computes the a Voronoi diagram on the unit sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHVORONOI.M - 10/10/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Convert to cartesian
    N_vert = size(dirs, 1);
    [tempx, tempy, tempz] = sph2cart(dirs(:,1), dirs(:,2), ones(N_vert,1));
    U_vert = [tempx, tempy, tempz];
    
    % number of triangles
    N_face = size(faces, 1);

    % Calculate the voronoi vertices for each triangle - for the unit sphere
    % are given by the unit normal vector of the triangle
    U_vor = zeros(N_face, 3);
    for n = 1:N_face
        i1 = faces(n, 1);
        i2 = faces(n, 2);
        i3 = faces(n, 3);

        r_12 = U_vert(i2, :) - U_vert(i1, :);
        r_13 = U_vert(i3, :) - U_vert(i1, :);
        r_normal = cross(r_12, r_13);

        u_normal = r_normal/norm(r_normal);

        U_vor(n, :) = u_normal;
    end
    
    voronoi.vert = U_vor;
    
    
    % Find duplicate vertices if any, due to two triangles sharing the same
    % circumscribed circle
    duplicates = zeros(N_face, 1);
    for n = 1:N_face
        if duplicates(n) == 0
            curVert = U_vor(n,:);

            for m = 1:N_face
                if n == m
                    m = m+1;
                elseif all(abs(curVert - U_vor(m,:))<1.0e-5)
                    duplicates(m) = n;
                else
                    m = m+1;
                end
            end
        end
    end
    
    
    % Calculate the voronoi polygons
    
    % find the an ordered sequence of the triangles around each vertex and
    % get the proper sequence of the voronoi vertices that constitute a
    % polygon
    for n = 1:N_vert
        faceIdx = []; % list of triangles that contain the specific vertex
        for m = 1:N_face
            if any(faces(m, :) == n)
                faceIdx(end+1) = m;
            end
        end
        
        % each triangle from the list contain the common vertex and two
        % other vertices - each triangle has two common vertices with each
        % other. One (brute) way of sorting the sequence is to pick one
        % triangle, find the neighbour triangle by finding their common
        % vertex, move to that triangle and iterate till all the number of
        % riangles have been checked.
        k = 1;
        currentfaceIdx = faceIdx(k); % pick-up one of the triangles in the list
        currentface = faces(currentfaceIdx, :); % the triangle's vertices
        currentvertIdx = find(currentface ~= n, 1); % pick-up one of the vertices that is not the central one
        currentvert = currentface(currentvertIdx);
        
        sorted = faceIdx(k); % this is the list that keeps the the ordered triangles
        
        notsorted = 1;
        
        while(notsorted)
        
            tempfacelist = faceIdx(find(currentfaceIdx ~= faceIdx)); % exclude the current triangle from the list

            for l = 1:length(tempfacelist)
                currentfaceIdx = tempfacelist(l);
                currentface = faces(currentfaceIdx, :);

                if any(currentface == currentvert) % if the currentvert exists in the current triangles vertices
                    sorted(end+1) = currentfaceIdx; % then it's the neighbour triangle - store its index
                    if length(sorted) == length(faceIdx) % if the sorted list has the length of faceIdx then done
                        notsorted = 0;
                        break
                    end

                    currentvertIdx = find((currentface ~= n).*(currentface ~= currentvert)); % find the next vertex from current triangle that excludes the central one and the pervious one
                    currentvert = currentface(currentvertIdx);
                    break
                end
            end
        
        end
        
        % remove the duplicate vertices from the list
        for i = 1:length(sorted)
            if duplicates(sorted(i)) ~= 0
                sorted(i) = duplicates(sorted(i));
            end
        end
        [~, notDupes] = unique(sorted, 'first');
        
        sorted = sorted(sort(notDupes));
        
        voronoi.face(n) = {sorted}; % save the sequence of voronoi vertices for each point
        sorted = [];
    end
    
end
