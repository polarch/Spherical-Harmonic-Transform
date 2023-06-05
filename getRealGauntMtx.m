function rGmtx = getRealGauntMtx(cGmtx)
%GETREALGAUNTMTX Construct a matrix of Gaunt coefficients for real SHs
%
% GETREALGAUNTMTX constructs the (N1+1)^2x(N2+1)^2x(N+1)^2 matrix of Gaunt 
% coefficients which represent the instegral of three real spherical  
% harmonics, such as
% G^q_{q',q''} = \int_\Omega Y_{q'}Y_{q''}Y^*_{q} \mathrm{d}\Omega.
%
% With Gaunt coefficients, the spherical harmonic coefficients of the
% product of two spherical functions can be given directly as a linear
% relationship between the harmonic coefficients of the two functions.
%
% The construction of those matrices is done through a transformation of 
% the respective matrices for complex SHs, since closed form solutions exist 
% only for the complex SH case.
%
%   cGmtx:  [(N1+1)^2 x (N2+1)^2 x (N+1)^2] tensor containing complex Gaunt  
% 	    coefficients with q' = n1^2+n1+m1+1 row indexing, q'' = n2^2+n2+m2+1
% 	    column indexing, and q = n^2+n+m+1 indexing at the third dimension.
%
%   rGmtx:  the output matrix of real Gaunt coefficients, with the same size 
% 	    and indexing as the input complex Gaunt tensor.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 5/6/2023
%   archontis.politis@atuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maximum SH orders of input complex Gaunt matrices
maxN = sqrt(size(cGmtx,3))-1;
maxN1 = sqrt(size(cGmtx,1))-1;
maxN2 = sqrt(size(cGmtx,2))-1;

% construct regular matrices transforming complex to real SHs
Wr2c1 = (complex2realSHMtx(N1)').';
Wr2c2 = (complex2realSHMtx(N2)').';
Wc2r = (complex2realSHMtx(N)').';

% transformation of the complex Gaunt matrices to real ones
rGmtx = zeros(size(cGmtx));
idx = 0;
for n=0:N
    W_n = Wc2r(idx+1:idx+2*n+1, idx+1:idx+2*n+1);
    for m=-n:n
        if m~=0
            w = W_n(m+n+1, [m+n+1 -m+n+1]);
            rGmtx(:,:,n^2+n+m+1) = Wr2c1.' * (w(1)*cGmtx(:,:,n^2+n+m+1) + w(2)*cGmtx(:,:,n^2+n-m+1)) * Wr2c2;
        elseif m==0
            w = W_n(n+1, n+1);
            rGmtx(:,:,n^2+n+m+1) = Wr2c1.' * cGmtx(:,:,n^2+n+m+1) * Wr2c2;
        end
    end
    idx = idx+2*n+1;
end

end
