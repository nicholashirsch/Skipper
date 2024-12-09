% ________________________________________________________________________
%
%                  STATE SPACE FORM OF EQUATIONS OF MOTION
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  D. DUENAS, N. HIRSCH, L. KATZ
% DATE:     12/03/2024
%
% OVERVIEW: THIS TAKES THE SOLVED EOMS AND CONSTRUCTS THE STATE SPACE MODEL
%           FOR THE SYSTEM. 

clear; clc;


% Define constants for the dynamic system as well as needed base points.
%            CM: Center of mass.
%            CT: Center of thrust.
%
%             M: Mass.
%             g: Acceleration due to gravity.
%          rho2: Distance between the CM and CT.
% Ixx, Iyy, Izz: Rotational inertia.
%            T0: Base point for thrust Taylor series approximation.
%         tauR0: Base point for RCS torque Taylor series approximation.
M = 1; g = 10; rho2 = 1;
Ixx = 1; Iyy = 1; Izz = 1;

T0    = 1;
tauR0 = 0;


A = [0 0 0 0 0            0          1 0 0 0 0 0 0;  %xDot
     0 0 0 0 0            0          0 1 0 0 0 0 0;  %yDot
     0 0 0 0 0            0          0 0 1 0 0 0 0;  %zDot
     0 0 0 0 0            0          0 0 0 1 0 0 0;  %phiDot 
     0 0 0 0 0            0          0 0 0 0 1 0 0;  %thetaDot
     0 0 0 0 0            0          0 0 0 0 0 1 0;  %psiDot
     0 0 0 0 T0/M         0          0 0 0 0 0 0 0;  %xDDot
     0 0 0 0 0            -T0/M      0 0 0 0 0 0 0;  %yDDot
     0 0 0 0 0            0          0 0 0 0 0 0 -g; %zDDot
     0 0 0 0 tauR0/Ixx    0          0 0 0 0 0 0 0;  %phiDDot
     0 0 0 0 0            -tauR0/Iyy 0 0 0 0 0 0 0;  %thetaDDot
     0 0 0 0 0            0          0 0 0 0 0 0 0;  %psiDDot
     0 0 0 0 0            0          0 0 0 0 0 0 0]; %0

B= [0   0    0     0; %xDot
    0   0    0     0; %yDot
    0   0    0     0; %zDot
    0   0    0     0; %phiDot
    0   0    0     0; %thetaDot
    0   0    0     0; %psiDot
    0   T0/M 0     0; %xDDot
    0   0    -T0/M 0; %yDDot
    1/M 0    0     0; %zDDot
    0   0    0     0; %phiDDot
    0   0    0     0; %thetaDDot
    0   0    0     0; %psiDDot
    0   0    0     1/Izz];%0

C = diag(ones(1, 13));

D = zeros(size(B));

states = {'x' 'y' 'z' 'phi' 'theta' 'psi' 'xDot' 'yDot' 'zDot' 'phiDot' 'thetaDot' 'psiDot' '1'};
inputs = {'T', 'xi', 'zeta', 'tauR'};
outputs = {'x' 'y' 'z' 'phi' 'theta' 'psi' 'xDot' 'yDot' 'zDot' 'phiDot' 'thetaDot' 'psiDot' '1'};

skipper = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs, 'name', 'Skipper')


