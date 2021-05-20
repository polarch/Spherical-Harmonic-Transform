function F = inverseSHT(F_N, dirs, basisType)
%INVERSE_SHT Perform the inverse spherical harmonic transform
%
%   F_N:    (N+1)^2 x L matrix of SH coefficients up to order N,  
%           with L spherical functions encoded as columns
%   dirs:   Kx2 matrix of directions that the inverse SHT is evaluated at, 
%           in [azimuth_1 inclination_1; ...; azimuth_K inclination_K] 
%           format (radians). Inclination is the polar angle from zenith
%           inclination = pi/2-elevation
%   basisType:  'complex' or 'real' spherical harmonics
%
%   F:      KxL matrix of reconstructed values at specified points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % vector of spherical harmonics
    N = sqrt(size(F_N,1)) - 1;
    Y_N = getSH(N, dirs, basisType);
    
    % perform the inverse transform up to degree N
    F = Y_N * F_N;

end
