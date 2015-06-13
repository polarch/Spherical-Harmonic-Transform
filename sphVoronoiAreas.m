function voronoi = sphVoronoiAreas(voronoi)
%SPHVORONOIAREAS Computes the areas of a voronoi diagram on the unit sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHVORONOIAREAS.M - 10/10/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vert = voronoi.vert;
    faces = voronoi.face;

    N_face = length(faces);

    areas = zeros(N_face, 1);
    for m = 1:N_face

        % current face
        face = faces{m};

        % number of vertices in the polygon
        N_poly = length(face);

        theta = zeros(N_poly, 1);
        for n = 1:N_poly
            % find vector between vertex origin and vertices 1 & 2
            r_01 = vert(face(1), :);
            r_02 = vert(face(2), :);

            % find normal vector to the great circle of 1 & 2
            r_2x1 = cross(r_02, r_01);
            
            % find tangent vector to great circle at 2
            r_21 = cross(r_2x1, r_02);

            % find vector between vertex origin and vertex 3
            r_03 = vert(face(3), :);

            % find normal vector to the great circle of 2 & 3
            r_2x3 = cross(r_02, r_03);
            
            % find tangent vector to great circle at 2
            r_23 = cross(r_2x3, r_02);            

            % normalise tangent vectors
            r_21 = r_21/norm(r_21);
            r_23 = r_23/norm(r_23);

            % get angle between the normals
            theta(n) = acos(dot(r_21, r_23));

            % shift the vertex list by one position and repeat
            face = circshift(face, [0 -1]);
        end

        % compute area of spherical polygon by spherical excess
        areas(m) = sum(theta) - (N_poly - 2)*pi;

    end
    
    voronoi.area = areas;

end

