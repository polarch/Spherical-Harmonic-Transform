function [F_N, Y_N] = directSHT(N, F, dirs, basisType, weights)
%LEASTSQUARESSHT Spherical harmonic transform of F for appropriate point sets
%
%   N:  maximum order of transform
%   F:  the spherical function evaluated at K directions 'dirs'
%       if F is a KxM matrix then the transform is applied to each column
%       of F, and returns the coefficients at each column of F_N
%       respectively
%   dirs:   [azimuth1 inclination1; ...; azimuthK inclinationK] angles in 
%           rads for each evaluation point, where inclination is the polar 
%           angle from zenith: inclination = pi/2-elevation
%   basisType:  'complex' or 'real' spherical harmonics
%   weights: vector of K weights for each direction, for weighted arrangements
%            such as Gauss-Legendre arrangements or Fliege-Maier sets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 7/2/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute the harmonic coefficients
    Y_N = getSH(N, dirs, basisType);
    Npoints = size(dirs,1);

    % non-weighted case for uniform arrangements
    if nargin<5 || ~exist('weights','var') || isempty(weights)
        
        % perform transform
        F_N = (4*pi/Npoints) * Y_N' * F;
    else
        
        % perform weights transform
        F_N = Y_N' * diag(weights) * F;
    end

end
