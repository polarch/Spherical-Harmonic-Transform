function x_rep = replicatePerOrder(x, dim)
%REPLICATEPERORDER Replicate l^th element 2*l+1 times across dimension
%
%   Replicates multidimensional array across dimension dim, so that the
%   l^th element of dim is replicated 2*l+1 times. that effectively has the
%   effect that the dimension grows from L to L^2 elements. This can be useful
%   in some spherical harmonic operations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 7/2/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    dim = 1;
end

ndimx = ndims(x);
sizx = size(x);
orderN = sizx(dim)-1;
sizx_rep = sizx;
sizx_rep(dim) = (orderN+1)^2;
x_rep = zeros(sizx_rep);

% I'm sure there's a better way to do this
str_dim_l = [];
rep_idx_l = [];
str_dim_r = [];
rep_idx_r = [];
for k=1:dim-1
    str_dim_l = [str_dim_l ':,'];
    rep_idx_l = [rep_idx_l '1,'];
end
for k=dim+1:ndimx
    str_dim_r = [str_dim_r ',:'];
    rep_idx_r = [rep_idx_r ',1'];
end
rep_idx = [rep_idx_l '2*n+1' rep_idx_r];

idx_n = 0;
for n=0:orderN
    eval(['x_n = x(' str_dim_l 'n+1' str_dim_r ');']);
    eval(['x_rep(' str_dim_l 'idx_n+(1:2*n+1)' str_dim_r ') = repmat(x_n,' rep_idx ');']);
    idx_n = idx_n+2*n+1;
end

end
