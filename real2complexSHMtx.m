function T_r2c = real2complexSHMtx(N)
%REAL2COMPLEXSHMTX Returns transformation matrix of real to complex SH
%   
%   Returns the unitary transformation matrix T_r2c the expresses the complex  
%   spherical harmonics with respect to the real ones, so that 
%   y_N = T_r2c * r_N, where r_N and y_N are the real and complex SH
%   vectors respectively, with all SHs up to order N included. For
%   normalisations and conventions used here for each base see the README
%   file.
%
%   N:      maximum order
%
%   T_r2c:  (N+1)^2 x (N+1)^2 basis transformation matrix
%           Note that this is not the matrix that converts the coefficients
%           That's (T_r2c').'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013, update 12/06/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maximum order
T_r2c = zeros((N+1)^2, (N+1)^2);
% n = 0
T_r2c(1,1) = 1;

idx = 1;
if N>0
    for n=1:N
        
        m = (1:n)';
        % form the diagonal
        diagT = [-1i*ones(n,1); sqrt(2)/2; (-1).^m]/sqrt(2);
        
        % form the antidiagonal
        adiagT = [ones(n,1); sqrt(2)/2; 1i*(-1).^m]/sqrt(2);
        
        % form the transformation matrix for the specific band n
        tempT = diag(diagT) + fliplr(diag(adiagT));
        
        T_r2c(idx + (1:2*n+1), idx + (1:2*n+1)) = tempT;
        idx = idx + 2*n+1;
    end
    
end

end
