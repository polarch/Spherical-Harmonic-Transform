function y_nm = sphConvolution(x_nm, h_n)
%SPHCONVOLUTION Performs spherical convolution of function X kernel H
%
%   The convolution of a spherical function x by a spherical filter h can
%   be performed very efficiently on the SH domain. From the filter kernel
%   h, only the m=0 degrees contribute to the output, and anyway it is
%   usually the case that the convolution kernel is axisymmetric, with
%   zero coefficients for m~=0.
%   
%   x_nm: (N+1)^2 coefficients of spherical function x(\theta,\phi)
%
%   h_n:  (N+1) coefficients of axisymmetric spherical kernel function 
%         h(\theta,\phi), with only SH coefficients for m=0 being non-zero.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% order
Nx = sqrt(length(x_nm)) - 1;
Nh = length(h_n) - 1;
N = min(Nx, Nh);

y_nm = zeros((N+1)^2, 1);
for n=0:N
    for m=-n:n
        q = n*(n+1)+m;
        y_nm(q+1) = 2*pi*sqrt(4*pi/(2*n+1))*x_nm(q+1)*h_n(n+1);
    end
end

end
