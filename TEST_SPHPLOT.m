%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 20/02/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test the three different functions for plotting 3-dimensional plots

% grid resolution
aziRes = 5;
polarRes = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CASE 1: Plot function from its coefficients

% plot a real 4th-order function
SHbasisType = 'complex';
F_N = getSH(4, [pi/4 pi/4], SHbasisType)';

h_ax = plotSphFunctionCoeffs(F_N, SHbasisType, aziRes, polarRes, 'real');

%% plot a complex 4-th order function (with random complex coefficients)
SHbasisType = 'complex';
F_N = randn(25,1) + 1i*randn(25,1);

h_ax = plotSphFunctionCoeffs(F_N, SHbasisType, aziRes, polarRes, 'complex');
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CASE 2: Plot function defined on a grid

% analytical expression for a real directional function 
% (4th-order cardioid)
fcosAlpha = @(azi, polar, azi0, polar0) cos(polar)*cos(polar0) + ...
    sin(polar)*sin(polar0).*cos(azi-azi0);  % function for dipole oriented at azi0, polar0
Nord = 4;
polar0 = pi/4;
azi0 = pi/4;
fcardioid = @(azi, polar) (1/2).^Nord * (1+fcosAlpha(azi, polar, azi0, polar0)).^Nord;

grid_dirs = grid2dirs(aziRes, polarRes);
F = fcardioid(grid_dirs(:,1), grid_dirs(:,2));
Fgrid = Fdirs2grid(F,aziRes,polarRes,1);

h_ax = plotSphFunctionGrid(Fgrid, aziRes, polarRes, 'real');

%% plot a complex 4-th order function defined on a grid 
% (e.g. a single complex spherical harmonic)
F = getSH(4, grid_dirs, 'complex');
F = F(:,20);
Fgrid = Fdirs2grid(F,aziRes,polarRes,1);

h_ax = plotSphFunctionGrid(Fgrid, aziRes, polarRes, 'complex');
axis equal
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CASE 3: Plot function defined on non-regular points

% get a unifrom arrangement of sampling points such as a tDesign
[~, tdesign_dirs] = getTdesign(21);
% convert from azi-elev to azi-inclination
tdesign_dirs = [tdesign_dirs(:,1) pi/2-tdesign_dirs(:,2)];

% analytical expression for a real directional function 
% (4rd-order dipole)
fcosAlpha = @(azi, polar, azi0, polar0) cos(polar)*cos(polar0) + ...
    sin(polar)*sin(polar0).*cos(azi-azi0);  % function for dipole oriented at azi0, polar0
Nord = 3;
polar0 = pi/4;
azi0 = pi/4;
fdipole = @(azi, polar) (fcosAlpha(azi, polar, azi0, polar0)).^Nord;

F = fdipole(tdesign_dirs(:,1), tdesign_dirs(:,2));

h_ax = plotSphFunctionTriangle(F, tdesign_dirs, 'real');

%% plot a complex 4-th order function defined on a grid 
% (just a single complex spherical harmonic)
F = getSH(4, tdesign_dirs, 'complex');
F = F(:,20);

h_ax = plotSphFunctionTriangle(F, tdesign_dirs, 'complex');
axis equal
colorbar
