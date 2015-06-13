function voronoi_weights = getVoronoiWeights(dirs)
%GETVORONOIWEIGHTS Compute voronoi diagram spherical areas for set of points
%
%   dirs:   directions of the K points on the sphere in the 
%           [azi1 elev1;... aziK elevK] convention
%
%   voronoi_weights:    Kx1 vector of the areas of the voronoi cells around
%                       each point. The areas sum to 4*pi.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 20/02/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dirs_x, dirs_y, dirs_z] = sph2cart(dirs(:,1), dirs(:,2), 1);
% perform delaunay triangulation
delaunay.vert = [dirs_x, dirs_y, dirs_z];
delaunay.face = sphDelaunay(dirs);
% get voronoi cells from delaunay
voronoi = sphVoronoi(dirs, delaunay.face);
% get areas of spherical voronoi polygons
voronoi = sphVoronoiAreas(voronoi);
voronoi_weights = voronoi.area;

end
