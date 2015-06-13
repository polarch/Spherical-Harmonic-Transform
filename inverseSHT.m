function F = inverseSHT(F_N, dirs, basisType)
%INVERSE_SHT Perform the inverse spherical harmonic transform
%
%   N:  maximum order of harmonics
%   F: the spherical function recreated at directions 'dirs'
%   dirs:   [azimuth inclination] angles in rads for each evaluation point,
%           where inclination is the polar angle from zenith
%           theta = pi/2-elevation
%   basisType:  'complex' or 'real' spherical harmonics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % vector of spherical harmonics
    N = sqrt(length(F_N)) - 1;
    Y_N = getSH(N, dirs, basisType);
    
    % perform the inverse transform up to degree N
    F = Y_N * F_N;

end
