% ________________________________________________________________________
%
%                        SKIPPER EQUATIONS OF MOTION
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  D. DUENAS
% DATE:     12/03/2024
%
% OVERVIEW: THIS TAKES THE SOLVED EOMs AND CONSTRUCTS THE STATE-SPACE MODEL
%           FOR THE SYSTEM
clear; clc;
import EOMAlgebra.*;
import EOMSims.*;

M = 1;
%zeta = 0;
%xi = 0;
%T = 0;
g = 10;
Ixx = 1;
Iyy = 1;
Izz = 1;
%tauR = 0;
T0 = 1;
tauR0 = 1;
rho2 = 1;


A = [0 1 0 0 0 0 0 0 0 0 0 0;   %*x = xDot
    0 0 0 0 0 0 0 0 g 0 0 0;    %*xDot = xDDot
    0 0 0 1 0 0 0 0 0 0 0 0;    %*y = yDot
    0 0 0 0 0 0 g 0 0 0 0 0;    %
    0 0 0 0 1 0 0 0 0 0 0 0;    %done
    0 0 0 0 0 0 0 0 0 0 0 0;    %done
    0 0 0 0 0 0 1 0 0 0 0 0;    %done
    0 0 0 0 0 0 0 0 0 0 0 0;    %done
    0 0 0 0 0 0 0 0 1 0 0 0;    %done
    0 0 0 0 0 0 0 0 0 0 0 0;    %done
    0 0 0 0 0 0 0 0 0 0 1 0;    %done
    0 0 0 0 0 0 0 0 0 0 0 0;];  %done

B= [0 0 0 0;                    %xDot
    0 T0/M 0 0;                   %xDDot
    0 0 0 0;                    %yDot
    0 0 -T0/M 0;                  %yDDot
    0 0 0 0;                    %zDot
    1/M 0 0 0;                  %zDDot how should I input gravity
    0 0 0 0;                    %phiDot
    0 0 0 0;                    %phiDDot
    0 0 0 0;                    %thetaDot
    rho2/Iyy 0 0 0;             %thetaDDot
    0 0 0 0;                    %psiDot
    0 0 T0*rho2/Izz 1/Izz;];    %psiDDot

C = [1 0 0 0 0 0 0 0 0 0 0 0;   %*x = xDot
    0 1 0 0 0 0 0 0 0 0 0 0;    %*xDot = xDDot
    0 0 1 0 0 0 0 0 0 0 0 0;    %*y = yDot
    0 0 0 1 0 0 0 0 0 0 0 0;    %done
    0 0 0 0 1 0 0 0 0 0 0 0;    %done
    0 0 0 0 0 1 0 0 0 0 0 0;    %done
    0 0 0 0 0 0 1 0 0 0 0 0;    %done
    0 0 0 0 0 0 0 1 0 0 0 0;    %done
    0 0 0 0 0 0 0 0 1 0 0 0;    %done
    0 0 0 0 0 0 0 0 0 1 0 0;    %done
    0 0 0 0 0 0 0 0 0 0 1 0;    %done
    0 0 0 0 0 0 0 0 0 0 0 1;];

D = [0 0 0 0;                    %xDot
    0 0 0 0;                   %xDDot
    0 0 0 0;                    %yDot
    0 0 0 0;                  %yDDot
    0 0 0 0;                    %zDot
    0 0 0 0;                  %zDDot how should I input gravity
    0 0 0 0;                    %phiDot
    0 0 0 0;                    %phiDDot
    0 0 0 0;                    %thetaDot
    0 0 0 0;             %thetaDDot
    0 0 0 0;                    %psiDot
    0 0 0 0;];    %psiDDot


states = {'x' 'xDot' 'y' 'yDot' 'z' 'zDot' 'phi' 'phiDot' 'theta' 'thetaDot' 'psi' 'psiDot'};
inputs = {'T' 'xi' 'zeta' 'tauR'};
outputs = {'x' 'xDot' 'y' 'yDot' 'z' 'zDot' 'phi' 'phiDot' 'theta' 'thetaDot' 'psi' 'psiDot'};
sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);