function R = euler2rotationMatrix(alpha, beta, gamma, convention)
%EULER2ROTATIONMATRIX Construct the rotation matrix from Euler angles
%
%   alpha:  first angle of rotation
%   beta:   second angle of rotation
%   gamma:  third angle of rotation
%
%   connvention: definition of the axes of rotation, e.g. for the y-convention
%                this should be 'zyz', for the x-convnention 'zxz', and for
%                the yaw-pitch-roll convention 'zyx'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 7/2/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rx = @(theta) [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
Ry = @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
Rz = @(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

switch convention(1)
    case 'x'
        R1 = Rx(alpha);
    case 'y'
        R1 = Ry(alpha);
    case 'z'
        R1 = Rz(alpha);
end

switch convention(2)
    case 'x'
        R2 = Rx(beta);
    case 'y'
        R2 = Ry(beta);
    case 'z'
        R2 = Rz(beta);
end

switch convention(3)
    case 'x'
        R3 = Rx(gamma);
    case 'y'
        R3 = Ry(gamma);
    case 'z'
        R3 = Rz(gamma);
end

R = R3*R2*R1;

end
