function w = w3j_stirling(j1, j2, j3, m1, m2, m3)
% W3J_STIRLING Computes Wigner-3j symbols using Stirling's approximation.
%
%   W3J_STIRLING computes the Wigner 3j symbol through the Racah formula found 
%   in   
%       http://mathworld.wolfram.com/Wigner3j-Symbol.html, Eq.7.
%
%   with the difference that it uses Stirling's approximation for large 
%   factorials ln(n!) = n*ln(n) - n, more on
%
%       http://mathworld.wolfram.com/StirlingsApproximation.html.
%
%   The code is based on the one by J.Pritchard found in:
%       http://massey.dur.ac.uk/jdp/code/w3j.m
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
    coeff = (-1)^(j1-j2-m3) * exp( 0.5*(logf(j1+m1) + logf(j1-m1) + logf(j2+m2) + ...
        logf(j2-m2) + logf(j3+m3) + logf(j3-m3) + logf(j1 + j2 - j3) + logf(j1 - j2 + j3) + ...
        logf(-j1 + j2 + j3) - logf(j1 + j2 + j3 + 1)));
    
    % summation over integers that do not result in negative factorials
    Sum_t = 0;
    for t = 0:N_t
        
        % check factorial for negative values, include in sum if not
        if j3-j2+t+m1 >= 0 && j3-j1+t-m2 >=0 && j1+j2-j3-t >= 0 && j1-t-m1 >=0 && j2-t+m2 >= 0
            x_t = exp( logf(t) + logf(j1+j2-j3-t) + logf(j3-j2+t+m1) + ...
                logf(j3-j1+t-m2) + logf(j1-t-m1) + logf(j2-t+m2) );
            Sum_t = Sum_t + (-1)^t/x_t;
        end
    end
    
    w = coeff*Sum_t;
end

end

% Stirling's large factorial approximation ln(n!) = nln(n)-n+0.5ln(2pin)
function y=logf(x)

if(x<170)
    y=log(factorial(x));
else
    y=x*log(x)-x+0.5*log(2*pi*x)+1/(12*x)-1/(360*x^3)+1/(1260*x^5)...
        -1/(1680*x^7)+1/(1188*x^9);
end

end