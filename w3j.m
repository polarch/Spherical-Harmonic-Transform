function w = w3j(j1, j2, j3, m1, m2, m3)
% W3J Computes Wigner-3j symbols.
%
%   W3J computes the Wigner 3j symbol through the Racah formula found in
%   
%   http://mathworld.wolfram.com/Wigner3j-Symbol.html, Eq.7.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 5/10/2014
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check selection rules (http://mathworld.wolfram.com/Wigner3j-Symbol.html)
if abs(m1)>abs(j1) || abs(m2)>abs(j2) || abs(m3)>abs(j3)
    w = 0;
elseif (m1+m2+m3~=0)
    w = 0;
elseif j3<abs(j1-j2) || j3>(j1+j2)  % triangle inequality
    w = 0;
    
else
    % evaluate the Wigner-3J symbol using the Racah formula (http://mathworld.wolfram.com/Wigner3j-Symbol.html)
    % number of terms for the summation
    N_t = max([j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3, j1+j2-j3, j2+j3-j1, j3+j1-j2]);
    
    % coefficients before the summation
    coeff1 = (-1)^(j1-j2-m3);
    coeff2 = factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)* ...
        factorial(j3+m3)*factorial(j3-m3);
    tri_coeff = factorial(j1 + j2 - j3)*factorial(j1 - j2 + j3)*factorial(-j1 + j2 + j3)/ ...
        factorial(j1 + j2 + j3 + 1);
    
    % summation over integers that do not result in negative factorials
    Sum_t = 0;
    for t = 0:N_t
        
        % check factorial for negative values, include in sum if not
        if j3-j2+t+m1 >= 0 && j3-j1+t-m2 >=0 && j1+j2-j3-t >= 0 && j1-t-m1 >=0 && j2-t+m2 >= 0
            x_t = factorial(t)*factorial(j1+j2-j3-t)*factorial(j3-j2+t+m1)* ...
                factorial(j3-j1+t-m2)*factorial(j1-t-m1)*factorial(j2-t+m2);
            Sum_t = Sum_t + (-1)^t/x_t;
        end
    end
    
    w = coeff1*sqrt(coeff2*tri_coeff)*Sum_t;
end

end

