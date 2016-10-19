function Wgrid = Fdirs2grid(W, aziRes, polarRes, CLOSED)
%FDIRS2GRID Replicate vector function values on a regular grid
%
%   Fdirs2grid takes a vector of values of function evaluated at a 
%   spherical grid with the grid2dirs function, and convert it back to a 
%   2D grid. Useful for plotting with functions such as surf, or for 
%   numerical integration numerically spherical functions.
%
%   W:  column vector of function values evaluated at each grid direction,
%       with the direction ordering given by grid2dirs. If W is a matrix
%       then each column is considered as a separate function to be
%       converted
%   aziRes: azimuth resolution of the grid in degrees (should be the same
%           as the one used in the grid2dirs function
%   polarRes:   inclination resolution of the grid in degrees (should be 
%               the same as the one used in the grid2dirs function
%   CLOSED: {0,1} if true then the returned matrix replicates the first
%           column of function values at 0deg azimuth also at 360deg,
%           useful for 3D plotting so that the shape does not have a
%           hole in the end (see plotSphFunction test script)
%
%   Wgrid:  if W is a vector then Wgrid is the 2D matrix of the function
%           values replicated on the grid points. If W is a matrix, then
%           Wgrid is a 3D matrix with one grid per column of W, on the 3rd
%           dimension.
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

if nargin<4
    CLOSED = 0;
end

Nphi = 360/aziRes;
Ntheta = 180/polarRes+1;

Nf = size(W, 2);
Wgrid = zeros(Nphi, Ntheta, Nf);
for i = 1:Nf
    
    Wgrid(:, 2:end-1, i) = reshape(W(2:end-1, i), Nphi, Ntheta-2);
    Wgrid(:, 1, i) = ones(Nphi, 1) * W(1, i);
    Wgrid(:, end, i) = ones(Nphi, 1) * W(end, i);
end

if Nf~=1
    Wgrid = permute(Wgrid, [2 1 3]);
else
    Wgrid = Wgrid.';
end

if CLOSED
    Wgrid = horzcat(Wgrid, Wgrid(:,1,:));
end
