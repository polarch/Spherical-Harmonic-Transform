function x_rep = replicatePerOrder(x)
%REPLICATEPERORDER Summary of this function goes here
%   Detailed explanation goes here

orderN = length(x)-1;

x_rep = [];
for n=0:orderN
    x_temp = ones(2*n+1,1)*x(n+1);
    x_rep = [x_rep; x_temp];
end


end

