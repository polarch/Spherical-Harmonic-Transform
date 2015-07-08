%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 20/02/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET-UP STUFF
%
% t-designs are uniform arrangements of points on the sphere that 
% fulfil exact integration of spherical polynnomials up to degree t, by
% simple summation of the values of the polynomial at these points. When 
% used for the spherical harmonic transform (SHT) up to order N, a design
% of N = floor(t/2) should be used, or equivalently t>=2N. For more details 
% see the function getTdesign().
%
% The Fliege-Maier nodes is another example of nearly-uniform arrangements
% that along with their respective integration weights can be used for
% direct integration through summation. When used for the SHT up to order
% N, the set with the index=N+1 should be loaded. For details see the 
% function getFliegeNodes().
%
% SHT with such special sampling arrangements, including other ones such as
% Gauss-Legendre quadratures, or Lebedev grids, then the function
% directSHT() can be used.
%
% When the function is sampled in non-uniform points, such as for example
% on a regular grid with equiangularly spaced points in azimuth and
% elevation which occurs often in practice, then it is better to perform 
% the SHT in the least-squares sense. Weights which can improve the
% conditioning can be passed to the function for a weighted least-squares
% solution, or they can be computed through the generic areas of voronoi
% cells around the sampling points, through the getVoronoiWeights()
% function.
%
% Examples of the above are presented below.


% Helper function to convert directions from Matlab's azimuth-elevation to
% azimuth-inclination
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

% Construct a band-limited spherical function for testing the SHT (4-th
% order cardioid function)
Nord = 4;
fcosAlpha = @(azi, polar, azi0, polar0) cos(polar)*cos(polar0) + ...
    sin(polar)*sin(polar0).*cos(azi-azi0);  % function for dipole oriented at azi0, polar0
% orientation of mainlobe
polar0 = pi/4;
azi0 = pi/4;
fcardioid = @(azi, polar) (1/2).^Nord * (1+fcosAlpha(azi, polar, azi0, polar0)).^Nord;

% evaluate function on a dense regular grid
dirs = grid2dirs(5, 5);
F = fcardioid(dirs(:,1), dirs(:,2));

% plot function using the spherical plotting function
figure
plotSphFunctionGrid(Fdirs2grid(F,5,5,1), 5,5,'real',gca);

% construct different sampling grids for the SHT of order 4
regular_grid = grid2dirs(30,30);
Kreg = size(regular_grid,1); % regular grid with 62 points

[~, tdesign_grid] = getTdesign(2*Nord);
tdesign_grid = aziElev2aziIncl(tdesign_grid); % convert to azi-incl
Ktdesign = size(tdesign_grid,1); % t-design of 36 points

[~, fliege_grid, fliege_weights] = getFliegeNodes(Nord+1);
fliege_grid = aziElev2aziIncl(fliege_grid); % convert to azi-incl
Kfliege = size(fliege_grid,1); % fliege nodes 25 points

% plot the different grids using the spherical plot functions
figure
h_ax = subplot(131); plotSphFunctionGrid(ones(length(0:30:180),length(0:30:360)), 30,30,'real', h_ax), title('Regular grid')
h_ax = subplot(132); plotSphFunctionTriangle(ones(length(tdesign_grid),1), tdesign_grid, 'real', h_ax), title('t-Design')
h_ax = subplot(133); plotSphFunctionTriangle(ones(length(fliege_grid),1), fliege_grid, 'real', h_ax), title('Maier-Fliege Design')

% get integration weights for the (non-uniform) regular grid, based on
% areas of voronoi cells around the sampling points
regular_weights = getVoronoiWeights(aziElev2aziIncl(regular_grid));

% check conditioning of the transform for up to order 5 (one more order than 
% the sampling schemes are for)
cond_N_reg = checkCondNumberSHT(Nord+1, regular_grid, 'complex', []);
cond_N_reg_weighted = checkCondNumberSHT(Nord+1, regular_grid, 'complex', regular_weights);
cond_N_tdes = checkCondNumberSHT(Nord+1, tdesign_grid, 'complex', []);
cond_N_fliege = checkCondNumberSHT(Nord+1, fliege_grid, 'complex', []);
cond_N_fliege_weighted = checkCondNumberSHT(Nord+1, fliege_grid, 'complex', fliege_weights);
figure
subplot(131), plot(0:Nord+1, [cond_N_reg cond_N_reg_weighted]), title('regular'), legend('unweighted','weighted')
subplot(132), plot(0:Nord+1, cond_N_tdes), title('t-Design')
subplot(133), plot(0:Nord+1, [cond_N_fliege cond_N_fliege_weighted]), title('Maier-Fliege design'), legend('unweighted','weighted')


%% PERFORM SHT

% sample function on the respective grids
F_reg = fcardioid(regular_grid(:,1), regular_grid(:,2));
F_tdes = fcardioid(tdesign_grid(:,1), tdesign_grid(:,2));
F_fliege = fcardioid(fliege_grid(:,1), fliege_grid(:,2));

% perform SHT for the regular grid using weighted least-squares
Fnm_reg = leastSquaresSHT(Nord, F_reg, regular_grid, 'complex', regular_weights);
% perform SHT for the t-Design with direct non-weighted summation
Fnm_tdes = directSHT(Nord, F_tdes, tdesign_grid, 'complex', []);
% perform SHT for the Fliege nodes using weighted summation
Fnm_fliege = directSHT(Nord, F_fliege, fliege_grid, 'complex', fliege_weights);

% plot coefficients
figure
plot(0:(Nord+1)^ 2-1, Fnm_reg)
hold on, plot(0:(Nord+1)^ 2-1, Fnm_tdes,'--r')
plot(0:(Nord+1)^ 2-1, Fnm_fliege, '-.*')
legend('regular','t-Design','Maier-Fliege')

% recreate the function at a very dense grid from the SH coefficients (using
% any of the above calculated coefficients) and plot
figure
plotSphFunctionCoeffs(Fnm_tdes, 'complex', 5, 5, 'real', gca)

% plot error between reconstructed/interpolated function at the grid points
% and the true values
Finterp = inverseSHT(Fnm_tdes, dirs, 'complex');
Ferror = abs(Finterp - F);
plotSphFunctionGrid(Fdirs2grid(Ferror,5,5,1), 5,5,'real'), title('Reconstruction error')


%% INTEGRATION OF SPHERICAL FUNCTIONS
%
% A test of the spectral theorem that states that the integral of the 
% product of two functions is equal to the dot product sum of their 
% respective spectral coefficients.

% Generate two random complex band-limited functions
Nf = 7;
Fnm = randn((Nf+1)^2,1) + 1i*randn((Nf+1)^2,1);
Ng = 5;
Gnm = randn((Ng+1)^2,1) + 1i*randn((Ng+1)^2,1);

% Generate a uniform grid for 12th-order integration
[~, tdesign_grid] = getTdesign(Nf+Ng);
Ktdes = size(tdesign_grid,1);

% Get the function values at the grid
Fgrid = inverseSHT(Fnm, aziElev2aziIncl(tdesign_grid), 'complex');
Ggrid = inverseSHT(Gnm, aziElev2aziIncl(tdesign_grid), 'complex');

% Integrate by direct summation
IntFG = (4*pi/Ktdes)*sum(Fgrid .* conj(Ggrid));

% Integrate directly using the SH coefficients
IntFG2 = Gnm(1:(min([Nf Ng])+1)^2)' * Fnm(1:(min([Nf Ng])+1)^2);

% Difference
disp(['Difference between integration in space and SH domain is ' num2str(IntFG - IntFG2)])


%% ROTATION OF AN AXISYMMETRIC SPHERICAL FUNCTION ON THE SH DOMAIN
%
% Rotation of general spherical functions directly by manipulation of their
% SH coefficients is fairly complicated, apart from the case of an
% axisymmetric function, which can be expressed as a weighted sum of
% Legendre polynomials (or a weighted sum of spherical harmonics of degree
% m=0). Such a function has N+1 non-zero SH coefficients. By rotating it to
% some arbitrary direction, all the SH coefficients are populated.

% define a real axisymmetric random pattern directly in the coefficient 
% domain (using only N+1 coefficients)
Nord = 6;
Fn = randn(Nord+1,1);

% target orientation to rotate the pattern
polar0 = pi/4;
azi0 = pi/4;
% fill the SH coefficients of the unrotated pattern with zeros for plotting
Fnm = rotateAxisCoeffs(Fn, 0,0, 'complex');
% SH coefficients of the rotated pattern
Fnm_rot = rotateAxisCoeffs(Fn, polar0, azi0, 'complex'); 

% plot the rotated and unrotated pattern
figure
h_ax = subplot(121); plotSphFunctionCoeffs(Fnm, 'complex', 5,5, 'real', h_ax), title('Non-rotated')
h_ax = subplot(122); plotSphFunctionCoeffs(Fnm_rot, 'complex', 5,5, 'real', h_ax), title('Rotated')


%% CONVOLUTION OF ONE FUNCTION WITH AN AXISYMMETRIC KERNEL
%
%   This script shows the result of the convolution of a spherical 
%   function x of order N=8 by a spherical filter h of lower order N=4. 
%   The original, filter and output functions are plotted. It's evident
%   that the result is of the lower order of the filter, which act as a
%   lowpass and eliminates the higher harmonics.

% generate a random complex 8th-order function
Nx = 8;
Xnm = randn((Nx+1)^2,1) + 1i*randn((Nx+1)^2,1);
% generate a random real 4th-order kernel
Nh = 4;
Hn = randn((Nh+1),1);
% fill the SH coefficients of the unrotated pattern with zeros for plotting
Hnm = rotateAxisCoeffs(Hn, 0,0, 'complex');

% perform convolution
Ynm = sphConvolution(Xnm, Hn);

% plot original function, kernel and output function
figure
h_ax = subplot(131); plotSphFunctionCoeffs(Xnm, 'complex', 5,5, 'complex', h_ax), title('Original')
h_ax = subplot(132); plotSphFunctionCoeffs(Hnm, 'complex', 5,5, 'real', h_ax), title('Kernel')
h_ax = subplot(133); plotSphFunctionCoeffs(Ynm, 'complex', 5,5, 'complex', h_ax), title('Convolution')


%% SOME MORE MANIPULATIONS OF SH COEFFICIENTS

%%% SH coefficients of conjugate function
%
%   Test that the SH coefficients of a function G, conjugate of F, are given 
%   correctly directly from the SH coefficients of F.

% generate a 4th-order complex function and get its values on a uniform
% grid
Nf = 6;
Fnm = randn((Nf+1)^2,1) + 1i*randn((Nf+1)^2,1);
[~, tdesign_grid] = getTdesign(2*Nf);
F_fliege = inverseSHT(Fnm, aziElev2aziIncl(tdesign_grid), 'complex');

% compute the SH coefficients of its conjugate function
G_fliege = conj(F_fliege);
Gnm = directSHT(Nf, G_fliege, aziElev2aziIncl(tdesign_grid), 'complex', []);

% compute coefficients of the conjugate function directly in the SH domain
Gnm_direct = conjCoeffs(Fnm);

% plot difference of the two
figure
plot(0:(Nf+1)^ 2-1, abs(Gnm - Gnm_direct))


%%% Conversion between real and complex SH coefficients
%
%   Test that the SH coefficients of a function, in the real and complex
%   SH basis, are given correctly by the complex-to-real and
%   real-to-complex transformation matrices.

% generate a random axisymmetric 4th-order function and rotate to some
% angle (polar0, azi0)
Nc = 4;
c_n = randn(Nc+1,1);
polar0 = pi/4;
azi0 = pi/4;

% get the SH coefficients of the rotated pattern into the two different
% bases
c_nm = rotateAxisCoeffs(c_n, polar0, azi0, 'complex');
r_nm = rotateAxisCoeffs(c_n, polar0, azi0, 'real');

% convert from one base to the other by conversion matrix
r_nm2 = complex2realCoeffs(c_nm);
c_nm2 = real2complexCoeffs(r_nm);

% plot difference
figure
plot(abs([(c_nm - c_nm2) (r_nm - r_nm2)])), legend('complex SH','real SH')


%% GAUNT COEFFICIENTS AND SPHERICAL MULTIPLICATION
%
% The Gaunt coefficient gives the integral of the product of three
% spherical harmonics and can be used in obtaining the coefficients of the
% product of two spherical functions directly from their respective
% coefficients

% Test that the coefficients of a product function of two spherical
% functions are equal to the ones derived directly by employing Gaunt
% coefficients (slow for high orders)

% Generate two random complex band-limited functions
Nf = 4;
Fnm = randn((Nf+1)^2,1) + 1i*randn((Nf+1)^2,1);
Ng = 2;
Gnm = randn((Ng+1)^2,1) + 1i*randn((Ng+1)^2,1);

% Generate a uniform grid for 12th-order integration
[~, tdesign_grid] = getTdesign(2*(Nf+Ng));

% Get the function values at the grid and that of the product function
Fgrid = inverseSHT(Fnm, aziElev2aziIncl(tdesign_grid), 'complex');
Ggrid = inverseSHT(Gnm, aziElev2aziIncl(tdesign_grid), 'complex');
Cgrid = Fgrid.*Ggrid;

% Get the SH coefficents of the product function by SHT
Cnm = directSHT(Nf+Ng, Cgrid, aziElev2aziIncl(tdesign_grid), 'complex', []);

% Get the coefficients of the product in the SH domain through Gaunt
% coefficients
Cnm_gaunt = sphMultiplication(Fnm, Gnm);

% Plot difference
plot(0:(Nf+Ng+1)^2-1, abs(Cnm - Cnm_gaunt))
