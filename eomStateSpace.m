% ________________________________________________________________________
%
%                  STATE SPACE FORM OF EQUATIONS OF MOTION
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  D. DUENAS, N. HIRSCH
% DATE:     12/03/2024
%
% OVERVIEW: THIS TAKES THE SOLVED EOMS AND CONSTRUCTS THE STATE SPACE MODEL
%           FOR THE SYSTEM. 


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
    
    B = zeros(12, 5);
    B(4, 1) = 1/M;
    B(4, 5) = g;
    B(5, 4) = T0/M;
    B(6, 3) = -T0/M;
    B(10, 2) = 1/Ixx;
    B(11, 3) = -T0*rho2/Iyy;
    B(12, 4) = -T0*rho2/Izz;
    
    C = diag(ones(1, 12));
    
    D = zeros(size(B));
    
    % Initialize LQR weight matricies.
    Q = diag(ones(1, 12));
    R = diag(ones(1, 5))
    R(5, 5) = 0.01; % Gravity isn't a real control.
    
    % Instantiate state space object. 
    sys = ss(A, B, C, D)
    
    % Solve for gain matrix, K, using 'lqr'.
    K = lqr(sys, Q, R)      




