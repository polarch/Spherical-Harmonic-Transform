function c_nm = sphMultiplication(a_nm, b_nm, G)
%SPHMULTIPLICATION Computes coefficients of a product spherical function
%   
%   The coefficients of a product function of two spherical functions can
%   be given directly as a linear combination of the coefficients of the
%   original functions, by employing Gaunt coefficients
%
%   a_nm: (N+1)^2 coefficients of spherical function a(\theta,\phi)
%   b_nm: (N+1)^2 coefficients of spherical function b(\theta,\phi)
%
%   c_nm:  (N+1)^2 coefficients of the product function 
%           c(\theta,\phi) = a(\theta,\phi)*b(\theta,\phi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% order
Na = sqrt(length(a_nm)) - 1;
Nb = sqrt(length(b_nm)) - 1;
Nc = Na+Nb;

c_nm = zeros((Nc+1)^2, 1);
% evaluate the coefficients of the product through the gaunt coefficients
if nargin<3
    G = gaunt_mtx(Na, Nb, Nc);
end
for n=0:Nc
    for m=-n:n
        q = n*(n+1)+m;
        c_nm(q+1) = a_nm.' * G(:,:,q+1) * b_nm;
        
    end
end

end
