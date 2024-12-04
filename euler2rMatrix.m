function rMatrix = euler2rMatrix(angle, axis)
% ________________________________________________________________________
%
%                            euler2rMatrix.m
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  N. HIRSCH
% DATE:     12/04/2024
%
% OVERVIEW: GENERATES A ROTATION MATRIX ABOUT THE GIVEN INTRINSIC AXIS BY THE
%           EULER GIVEN ANGLE.
% INPUTS:  
%   angle: Euler angle of rotation - double or symbolic.
%    axis: Intrinsic axis to peform rotation about - 3, 2, or 1.
% OUTPUTS:
%   rMatrix: Rotation matrix.

switch axis
    case 1 % x-axis rotation.
        rMatrix = [1       0           0    ;
                   0  cos(angle)  sin(angle);
                   0 -sin(angle)  cos(angle)];
    case 2 % y-axis rotation.
         rMatrix = [cos(angle)      0    -sin(angle);
                        0           1          0    ;
                    sin(angle)      0     cos(angle)];
    case 3 % z-axis rotation.
         rMatrix = [cos(angle)  sin(angle) 0;
                   -sin(angle)  cos(angle) 0;
                        0           0     1];
end