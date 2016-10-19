function dirs = grid2dirs(aziRes, polarRes, POLAR_OR_ELEV, ZEROED_OR_CENTERED)
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
    error('azimuth or elevation resolution should divide exactly 360 and 180deg')
end

switch nargin
    case 3
        ZEROED_OR_CENTERED = 1;
    case 2
        ZEROED_OR_CENTERED = 1;        
        POLAR_OR_ELEV=1;
    case {1,0}
        error('Not enough arguments');
end

if ZEROED_OR_CENTERED, phi = (0:aziRes:360-aziRes)*pi/180;
else phi = (-180:aziRes:180-aziRes)*pi/180;
end

if POLAR_OR_ELEV, theta = (0:polarRes:180)*pi/180; % polar angle
else theta = (-90:polarRes:90)*pi/180; % elevation angle
end
Nphi = length(phi);
Ntheta = length(theta);

for i=2:Ntheta-1   
    dirs((i-2)*Nphi + (1:Nphi) + 1, :) = [phi' theta(i)*ones(size(phi'))];
end

if POLAR_OR_ELEV
    dirs(1, :) = [0 0];
    dirs(end+1,:) = [0 pi];    
else
    dirs(1, :) = [0 -pi/2];
    dirs(end+1,:) = [0 pi/2];    
end
