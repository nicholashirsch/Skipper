function [A, B, C, D, K] = lqrMatrices(M, g, rho2, Ixx, Iyy, Izz, T0, tauR0)
% ________________________________________________________________________
%
%                              LQR MATRICES
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  N. HIRSCH
% DATE:     1/20/25
%
% OVERVIEW: Calculates A, B, C, D, Q, and R matices for use in a state
% space Simulink model and LQR controller.

    clc;

    % Initialize state space matricies. Since they are sparse simply use zeros
    % and then input non-zero terms individually.
    A = zeros(12, 12);
    A(1, 4) = 1;
    A(2, 5) = 1;
    A(3, 6) = 1;
    A(5, 9) = g;
    A(6, 8) = g;
    A(7, 10) = 1;
    A(8, 11) = 1;
    A(9, 12) = 1;
    
    B = zeros(12, 4);
    B(4, 1) = 1/M;
    B(5, 4) = T0/M;
    B(6, 3) = -T0/M;
    B(10, 2) = 1/Ixx;
    B(11, 3) = -T0*rho2/Iyy;
    B(12, 4) = -T0*rho2/Izz;
    
    C = diag(ones(1, 12));
    
    D = zeros(size(B));
    
    % Initialize LQR weight matricies.
    Q = diag(ones(1, 12));
    R = diag(ones(1, 4));
    
    % Instantiate state space object. 
    sys = ss(A, B, C, D);
    
    % Solve for gain matrix, K, using 'lqr'.
    K = lqr(sys, Q, R);
end
