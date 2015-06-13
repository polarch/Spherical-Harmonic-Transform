function faces = sphDelaunay(dirs)
%SPHDELAUNAY Computes the Delaunay triangulation on the unit sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHDELAUNAY.M - 10/10/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to cartesian
N_vert = size(dirs, 1);
[tempx, tempy, tempz] = sph2cart(dirs(:,1), dirs(:,2), ones(N_vert,1));
U_vert = [tempx, tempy, tempz];

% Find the convex hull of the points on the sphere - in this special case
% the result equals the Delaunay triangulation of the points
faces = convhulln(U_vert);

% Invert the triangles
faces = faces(:, 3:-1:1);

% Shift the results to begin each triangle from the smallest entry
for n = 1:size(faces,1)
    tempface = faces(n,:);
    [~, minIdx] = min(tempface);
    faces(n, :) = circshift(tempface, [0 1-minIdx]);
end

% Sort through triangles with smaller entries first
faces = sortrows(faces, 1); % sort through first entry
maxentry = max(faces(:,1)); % sort through second entry
n = 1;
while n <= maxentry
    startIdx = find(faces(:,1) == n, 1, 'first');
    if ~isempty(startIdx)
        endIdx = find(faces(:,1) == n, 1, 'last');
        faces(startIdx:endIdx, :) = sortrows(faces(startIdx:endIdx, :), 2);
        n = n + 1;
    else
        n = n + 1;
    end
end
