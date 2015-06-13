function g_nm = conjCoeffs(f_nm)
%CONJCOEFFS Get the complex SH coefficients of a conjugate sph. function
%
%   f_nm:   (N+1)^2 coefficients of original spherical function f(\theta,\phi)
%
%   g_nm:   (N+1)^2 coefficients of conjugate spherical function 
%               g(\theta,\phi) = (f(\theta,\phi))*
%
%   The conversion is based on the relation Y^*_{nm} = (-1)^m Y_{n(-m)}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% order
N = sqrt(length(f_nm)) - 1;

g_nm = zeros((N+1)^2, 1);
for n=0:N
    for m=-n:n
        q_g = n*(n+1)+m;
        q_f = n*(n+1)-m;
        g_nm(q_g+1) = (-1)^m * conj(f_nm(q_f+1));
    end
end


end
