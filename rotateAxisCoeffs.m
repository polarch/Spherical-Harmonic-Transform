function c_nm = rotateAxisCoeffs(c_n, theta_0, phi_0, basisType)
%ROTATEAXISCOEFFS Get spherical coefficients for a rotated axisymmetric pattern
%
%   c_n: N+1 coefficients describing a rotationally symmetric pattern of
%        order N, expressed as a sum of spherical harmonics of degree m=0
%        (sum of Legendre polynomials)
%   theta_0: polar rotation for the pattern
%   phi_0: azimuthal rotation for the pattern
%   basisType:  complex or real spherical harmonics
%
%   c_nm: (N+1)^2 coefficients of rotated pattern expressed as a sum of
%         spherical harmonics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(c_n)-1;

c_nm = zeros((N+1)^2, 1);
Y_N = getSH(N, [phi_0 theta_0], 'complex');
for n=0:N
    m = -n:n;
    q = n*(n+1)+m;
    c_nm(q+1) = sqrt(4*pi/(2*n+1)) * c_n(n+1) * conj(Y_N(q+1));
end

% convert to real SH coefficients if asked
if isequal(basisType, 'real'), c_nm = complex2realCoeffs(c_nm);
end
