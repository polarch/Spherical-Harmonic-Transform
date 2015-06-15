%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/06/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Test for rotation of spherical functions directly on the spherical
%   harmonic domain.
%

% resolution of the regular grid to plot the function
aziRes = 10;
polarRes = 10;

% desired orientation of the rotated function
alpha = pi/2;
beta = pi/2;
gamma = pi/4;

% define the rotation matrices succesively for plotting
R1 = euler2rotationMatrix(alpha, 0, 0, 'zyz');
R2 = euler2rotationMatrix(alpha, beta, 0, 'zyz');
R3 = euler2rotationMatrix(alpha, beta, gamma, 'zyz');

% generate random 3rd-order real function (directly in the SHD)
N = 3;
B_N = randn((N+1)^2,1);

% generate random 3rd-order complex function (directly in the SHD)
C_N = randn((N+1)^2,1) + 1i*randn((N+1)^2,1);

%% Rotate the function in the space domain directly as reference for comparison

aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; % convert elevation to inclinations due to definition of spherical harmonics

% rotate the directions of the evaluation grid according to the euler
% angles
dirs = grid2dirs(aziRes, polarRes);
[U_dirs(:,1),U_dirs(:,2),U_dirs(:,3)] = sph2cart(dirs(:,1), pi/2-dirs(:,2), 1);
U_rot1 = U_dirs * R1.';
[dirs1(:,1),dirs1(:,2)] = cart2sph(U_rot1(:,1), U_rot1(:,2), U_rot1(:,3));
U_rot2 = U_dirs * R2.';
[dirs2(:,1),dirs2(:,2)] = cart2sph(U_rot2(:,1), U_rot2(:,2), U_rot2(:,3));
U_rot3 = U_dirs * R3.';
[dirs3(:,1),dirs3(:,2)] = cart2sph(U_rot3(:,1), U_rot3(:,2), U_rot3(:,3));
% compute the functions at the original directions 
B0 = inverseSHT(B_N, dirs, 'real');
C0 = inverseSHT(C_N, dirs, 'complex');

% plot the real function before and after rotation
h1 = figure;
subplot(241)
plotSphFunctionCoeffs(B_N,'real',5,5,'real',gca)
subplot(242)
plotSphFunctionTriangle(B0,aziElev2aziIncl(dirs1),'real',gca)
subplot(243)
plotSphFunctionTriangle(B0,aziElev2aziIncl(dirs2),'real',gca)
subplot(244)
plotSphFunctionTriangle(B0,aziElev2aziIncl(dirs3),'real',gca)

% plot the complex function before and after rotation
h2 = figure;
subplot(241)
plotSphFunctionCoeffs(C_N,'complex',5,5,'complex',gca)
subplot(242)
plotSphFunctionTriangle(C0,aziElev2aziIncl(dirs1),'complex',gca)
subplot(243)
plotSphFunctionTriangle(C0,aziElev2aziIncl(dirs2),'complex',gca)
subplot(244)
plotSphFunctionTriangle(C0,aziElev2aziIncl(dirs3),'complex',gca)

%%

% rotate real coefficients
R_rSH1 = getSHrotMtx(R1, N, 'real');
R_rSH2 = getSHrotMtx(R2, N, 'real');
R_rSH3 = getSHrotMtx(R3, N, 'real');

% plot before and after rotation
figure(h1)
subplot(245)
plotSphFunctionCoeffs(B_N,'real',5,5,'real',gca)
subplot(246)
plotSphFunctionCoeffs(R_rSH1*B_N ,'real',5,5,'real',gca)
subplot(247)
plotSphFunctionCoeffs(R_rSH2*B_N,'real',5,5,'real',gca)
subplot(248)
plotSphFunctionCoeffs(R_rSH3*B_N,'real',5,5,'real',gca)

%%

% rotate complex coefficients
D_cSH1 = getSHrotMtx(R1, N, 'complex');
D_cSH2 = getSHrotMtx(R2, N, 'complex');
D_cSH3 = getSHrotMtx(R3, N, 'complex');

% plot before and after rotation
figure(h2)
subplot(245)
plotSphFunctionCoeffs(C_N,'complex',5,5,'complex',gca)
subplot(246)
plotSphFunctionCoeffs(D_cSH1*C_N,'complex',5,5,'complex',gca)
subplot(247)
plotSphFunctionCoeffs(D_cSH2*C_N,'complex',5,5,'complex',gca)
subplot(248)
plotSphFunctionCoeffs(D_cSH3*C_N,'complex',5,5,'complex',gca)
