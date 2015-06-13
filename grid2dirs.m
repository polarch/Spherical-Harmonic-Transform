function dirs = grid2dirs(aziRes, polarRes)
%GRID2DIRS Create a vector of spherical grid points
%
%   Grid2dirs creates a vector of a spherical grid directions, based on a 
%   given azimuthal and polar resolution. Useful for evaluating a spherical
%   function on a regular grid of points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(360, aziRes) ~= 0 || mod(180, polarRes) ~= 0
    error('azimuth or elevation resolution should divide exactly 180 and 150deg')
end

phi = (0:aziRes:360-aziRes)*pi/180;
theta = (0:polarRes:180)*pi/180;
Nphi = length(phi);
Ntheta = length(theta);

dirs(1, :) = [0 0];
for i=2:Ntheta-1
    
    dirs((i-2)*Nphi + (1:Nphi) + 1, :) = [phi' theta(i)*ones(size(phi'))];
end

dirs(end+1,:) = [0 pi];
