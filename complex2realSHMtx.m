function T_c2r = complex2realSHMtx(N)
%COMPLEX2REALSHMTX Returns transformation matrix of complex to rel SH
%   
%   Returns the unitary transformation matrix T_c2r the expresses the real  
%   spherical harmonics with respect to the complex ones, so that 
%   r_N = T_c2r * y_N, where r_N and y_N is are the real and complex SH
%   vectors respectively, with all SHs up to order N included. For
%   normalisations and conventions used here for each base see the README
%   file.
%
%   N:      maximum order
%
%   T_c2r:  (N+1)^2 x (N+1)^2 basis transformation matrix
%           Note that this is not the matrix that converts the coefficients
%           That's (T_c2r').'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013, update 12/06/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maximum order
T_c2r = zeros((N+1)^2, (N+1)^2);
% n = 0
T_c2r(1,1) = 1;

idx = 1;
if N>0
    for n=1:N
        
        m = (1:n)';
        % form the diagonal
        diagT = [1i*ones(n,1); sqrt(2)/2; (-1).^m]/sqrt(2);
        
        % form the antidiagonal
        adiagT = [-1i*(-1).^m(end:-1:1); sqrt(2)/2; ones(n,1)]/sqrt(2);
        
        % form the transformation matrix for the specific band n
        tempT = diag(diagT) + fliplr(diag(adiagT));
        
        T_c2r(idx + (1:2*n+1), idx + (1:2*n+1)) = tempT;
        idx = idx + 2*n+1;
    end
    
end

end
