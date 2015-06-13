function Lx = legendre2(N, x)
%LEGENDRE2 Same as matlab Legendre extended for negative orders m<0.
%
%   Pn(-m)(x) = (-1)^m * (n-m)!/(n+m)! * Pnm(x), for -m<0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEGENDRE2.M - 10/10/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = legendre(N, x);

if N ~= 0
    
    m = (0:N)';
    norm = (-1).^m .* factorial(N-m)./factorial(N+m);

    if isrow(x)
        dimx = fliplr(size(x));
    else
        dimx = size(x);
    end
    Lx_neg = repmat(norm, [1 dimx]) .* Lx;

    % I'm sure there a better way for this
    ndimx = ndims(x);
    str_dims = [];
    for k=1:ndimx
        str_dims = [str_dims ', :'];
    end
    % throw away first row and flip row order
    eval(['Lx_neg = Lx_neg(end:-1:2' str_dims ');'])
    
    Lx = [Lx_neg; Lx];
end
