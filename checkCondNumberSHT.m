function cond_N = checkCondNumberSHT(N, dirs, basisType, W)
%CHECKCONDNUMBER  Computes the condition number for a least-squares SHT
%
%   N:  maximum order to be tested for the given set of points
%   dirs:   [azimuth inclination] angles in rads for each evaluation point,
%           where inclination is the polar angle from zenith
%           theta = pi/2-elevation
%   basisType:  'complex' or 'real' spherical harmonics
%   W:  weights for each measurement point to condition the inversion, in
%       case of a weighted least-squares (Ndirs x 1 vector)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute the harmonic coefficients
    Y_N = getSH(N, dirs, basisType);
    
    % compute condition number for progressively increasing order up to N
    cond_N = zeros(N+1,1);
    for n = 0:N
        Y_n = Y_N(:, 1:(n+1)^2);
        if ~exist('W','var') || isempty(W)
            YY_n = (Y_n' * Y_n);
        else
            YY_n = (Y_n' * diag(W) * Y_n);
        end
        cond_N(n+1) = cond(YY_n);
    end

end
