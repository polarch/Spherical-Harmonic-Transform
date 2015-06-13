function C_N = real2complexCoeffs2(R_N)
%REAL2COMPLEXCOEFFS Convert SH coeffs from the real to complex basis
%   
%   Converts the vector of (N+1)^2xK SH coefficients of K functions on the
%   real SH base, to the respective ones of the complex SH base. For
%   normalisations and conventions used here for each base see the README
%   file.
%
%   R_N:  matrix of (N+1)^2 x K real SH coefficients
%
%   C_N:  matrix of (N+1)^2 x K complex SH coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013, update 12/06/2014
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maximum order
N = sqrt(size(R_N, 1)) -1;
T_c2r = complex2realCoeffsMtx2(N);

% get real-to-complex basis tranformation matrix by hermitian transpose of
% the unitary complex-to-real matrix computed above
T_r2c = T_c2r';
% convert coefficients
C_N = (T_r2c').' * R_N;
    
end
