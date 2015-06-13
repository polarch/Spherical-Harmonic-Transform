function Y_N = getSH(N, dirs, basisType)
%GETSH Get Spherical harmonics up to order N
%
%   N:  maximum order of harmonics
%   dirs:   [azimuth_1 inclination_1; ...; azimuth_K inclination_K] angles 
%           in rads for each evaluation point, where inclination is the 
%           polar angle from zenith: inclination = pi/2-elevation
%   basisType:  'complex' or 'real' spherical harmonics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ndirs = size(dirs, 1);
    Nharm = (N+1)^2;

    Y_N = zeros(Nharm, Ndirs);
    idx_Y = 0;
    for n=0:N
        
        m = (0:n)';
        if isequal(basisType, 'complex')
            % vector of unnormalised associated Legendre functions of current order
            Lnm = legendre(n, cos(dirs(:,2)'));
            
            % normalisations
            norm = sqrt( (2*n+1)*factorial(n-m) ./ (4*pi*factorial(n+m)) );
            
            % convert to matrix, for direct matrix multiplication with the rest
            Nnm = norm * ones(1,Ndirs);
            
            % spherical harmonics of current order
            Exp = exp(1i*m*dirs(:,1)');
            Ynm_pos = Nnm .* Lnm .* Exp;
            
            if n~=0
                % conjugate relation Y_{n(-m)} = (-1)^m (Y_{nm})^* for negative
                % degrees
                condon = (-1).^m(end:-1:2) * ones(1,Ndirs);
                Ynm_neg = condon .* conj(Ynm_pos(end:-1:2,:));
            else
                Ynm_neg = [];
            end
            
            % final SHs
            Ynm = [Ynm_neg; Ynm_pos];
            
        elseif isequal(basisType, 'real')
            % vector of unnormalised associated Legendre functions of current order
            Lnm_real = legendre(n, cos(dirs(:,2)'));
            if n~=0
                % cancel the Condon-Shortley phase from the definition of
                % the Legendre functions to result in signless real SH
                condon = (-1).^[m(end:-1:2);m] * ones(1,Ndirs);
                Lnm_real = condon .* [Lnm_real(end:-1:2, :); Lnm_real];
            end
            
            % normalisations
            norm_real = sqrt( (2*n+1)*factorial(n-m) ./ (4*pi*factorial(n+m)) );
            
            % convert to matrix, for direct matrix multiplication with the rest
            Nnm_real = norm_real * ones(1,Ndirs);
            if n~=0
                Nnm_real = [Nnm_real(end:-1:2, :); Nnm_real];
            end            
            
            CosSin = zeros(2*n+1,Ndirs);
            % zero degree
            CosSin(n+1,:) = ones(1,size(dirs,1));
            % positive and negative degrees
            if n~=0
                CosSin(m(2:end)+n+1,:) = sqrt(2)*cos(m(2:end)*dirs(:,1)');
                CosSin(-m(end:-1:2)+n+1,:) = sqrt(2)*sin(m(end:-1:2)*dirs(:,1)');
            end
            Ynm = Nnm_real .* Lnm_real .* CosSin;
            
        end
        
        Y_N(idx_Y+1:idx_Y+(2*n+1), :) = Ynm;
        idx_Y = idx_Y + 2*n+1;
    end
    
    % transpose to Ndirs x Nharm
    Y_N = Y_N.';
    
end
