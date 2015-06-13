function [F_N, Y_N] = leastSquaresSHT(N, F, dirs, basisType, weights)
%LEASTSQUARESSHT Spherical harmonic transform of F using least-squares
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
%   weights:    vector of K weights for each direction, for weighted
%               least-squares
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013, updated 7/2/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % compute the harmonic coefficients
    Y_N = getSH(N, dirs, basisType);
    Npoints = size(dirs,1);
    
    % non-weighted case for uniform arrangements
    if nargin<5 || ~exist('weights','var') || isempty(weights)
        
        % perform transform in the least squares sense
        F_N = pinv(Y_N)*F;
    else
        
        % perform weighted least-squares transform
        F_N = (Y_N'*diag(weights)*Y_N) \ (Y_N'*diag(weights) * F);
    end        
        

end
