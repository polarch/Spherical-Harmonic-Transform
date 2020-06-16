function R = getSHrotMtx(Rxyz, L, basisType)
%GETSHROTMTX Rotation matrices of real/complex spherical harmonics
%   GETSHROTMTX computes the rotation matrices for real spherical
%   harmonics according to the recursive method of Ivanic and Ruedenberg,
%   as can be found in
%
%       Ivanic, J., Ruedenberg, K. (1996). Rotation Matrices for Real 
%       Spherical Harmonics. Direct Determination by Recursion. The Journal 
%       of Physical Chemistry, 100(15), 6342?6347.
%
%   and with the corrections of
%
%       Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real 
%       Spherical Harmonics. Direct Determination by Recursion Page: Additions 
%       and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100.
%
%   The code implements directly the equations of the above publication and
%   is based on the code, with modifications, found in submision
%   http://www.mathworks.com/matlabcentral/fileexchange/15377-real-valued-spherical-harmonics
%   by Bing Jian
%   and the C++ implementation found at
%   http://mathlib.zfx.info/html_eng/SHRotate_8h-source.html.
%
%   Apart from the real rotation matrices, for which the algorithm is
%   defined, the function returns also the complex SH rotation matrices, by
%   using the transformation matrices from real to complex spherical
%   harmonics. This way bypasses computation of the complex rotation
%   matrices based on explicit formulas of Wigner-D matrices, which would
%   have been much slower and less numerically robust.
%
%   TODO: Recursive algorithm directly for complex SH as found in 
%       
%       Choi, C. H., Ivanic, J., Gordon, M. S., & Ruedenberg, K. (1999). Rapid 
%       and stable determination of rotation matrices between spherical 
%       harmonics by direct recursion. The Journal of Chemical Physics, 
%       111(19), 8825.
%
%   Inputs:
%
%       L: the maximum band L up to which the rotation matrix is computed
%       Rxyz: the normal 3x3 rotation matrix for the cartesian system
%       basisType: 'real' or 'complex' SH
%
%   Outputs:
%
%       R: the (L+1)^2x(L+1)^2 block diagonal matrix that rotates the frame
%       to the desired orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/06/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    basisType = 'real';
end
% allocate total rotation matrix
R = zeros((L+1)^2);

% initialize zeroth and first band rotation matrices for recursion
%	Rxyz = [Rxx Rxy Rxz
%           Ryx Ryy Ryz
%           Rzx Rzy Rzz]
%
% zeroth-band (l=0) is invariant to rotation
R(1) = 1;

% the first band (l=1) is directly related to the rotation matrix
R_1(-1+2,-1+2) = Rxyz(2,2);
R_1(-1+2, 0+2) = Rxyz(2,3);
R_1(-1+2, 1+2) = Rxyz(2,1);
R_1( 0+2,-1+2) = Rxyz(3,2);
R_1( 0+2, 0+2) = Rxyz(3,3);
R_1( 0+2, 1+2) = Rxyz(3,1);
R_1( 1+2,-1+2) = Rxyz(1,2);
R_1( 1+2, 0+2) = Rxyz(1,3);
R_1( 1+2, 1+2) = Rxyz(1,1);

R(2:4,2:4) = R_1;
R_lm1 = R_1;

% compute rotation matrix of each subsequent band recursively
band_idx = 4;
for l = 2:L
    
    R_l = zeros(2*l+1);    
    for m=-l:l
        for n=-l:l
            
            % compute u,v,w terms of Eq.8.1 (Table I)
            d = (m==0); % the delta function d_m0
            if abs(n)==l
                denom = (2*l)*(2*l-1);
            else
                denom = (l*l-n*n);
            end
            u = sqrt((l*l-m*m)/denom);
            v = sqrt((1+d)*(l+abs(m)-1)*(l+abs(m))/denom)*(1-2*d)*0.5;
            w = sqrt((l-abs(m)-1)*(l-abs(m))/denom)*(1-d)*(-0.5);
            
            % computes Eq.8.1
            if u~=0, u = u*U(l,m,n,R_1,R_lm1); end
            if v~=0, v = v*V(l,m,n,R_1,R_lm1); end
            if w~=0, w = w*Wf(l,m,n,R_1,R_lm1); end
            R_l(m+l+1,n+l+1) = u + v + w;
        end
    end
    R(band_idx+(1:2*l+1), band_idx+(1:2*l+1)) = R_l;
    R_lm1 = R_l;
    band_idx = band_idx + 2*l+1;
end

% if the rotation matrix is needed for complex SH, then get it from the one
% for real SH by the real-to-complex-transformation matrices
if isequal(basisType, 'complex')
    W = complex2realSHMtx(L);
    R = W.'*R*conj(W);
end

end


% functions to compute terms U, V, W of Eq.8.1 (Table II)
function [ret] = U(l,m,n,R_1,R_lm1)

ret = P(0,l,m,n,R_1,R_lm1);

end

function [ret] = V(l,m,n,R_1,R_lm1)

if (m==0)
    p0 = P(1,l,1,n,R_1,R_lm1);
    p1 = P(-1,l,-1,n,R_1,R_lm1);
    ret = p0+p1;
else
    if (m>0)
        d = (m==1);
        p0 = P(1,l,m-1,n,R_1,R_lm1);
        p1 = P(-1,l,-m+1,n,R_1,R_lm1);        
        ret = p0*sqrt(1+d) - p1*(1-d);
    else
        d = (m==-1);
        p0 = P(1,l,m+1,n,R_1,R_lm1);
        p1 = P(-1,l,-m-1,n,R_1,R_lm1);        
        ret = p0*(1-d) + p1*sqrt(1+d);
    end
end

end

function [ret] = Wf(l,m,n,R_1,R_lm1)

if (m==0)
    error('should not be called')
else
    if (m>0)
        p0 = P(1,l,m+1,n,R_1,R_lm1);
        p1 = P(-1,l,-m-1,n,R_1,R_lm1);        
        ret = p0 + p1;
    else
        p0 = P(1,l,m-1,n,R_1,R_lm1);
        p1 = P(-1,l,-m+1,n,R_1,R_lm1);        
        ret = p0 - p1;
    end
end

end

% function to compute term P of U,V,W (Table II)
function [ret] = P(i,l,a,b,R_1,R_lm1)

ri1 = R_1(i+2,1+2);
rim1 = R_1(i+2,-1+2);
ri0 = R_1(i+2,0+2);

if (b==-l)
    ret = ri1*R_lm1(a+l,1) + rim1*R_lm1(a+l, 2*l-1);
else
    if (b==l)
        ret = ri1*R_lm1(a+l,2*l-1) - rim1*R_lm1(a+l, 1);        
    else
        ret = ri0*R_lm1(a+l,b+l);
    end
end

end
