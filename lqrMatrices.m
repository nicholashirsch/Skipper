function [A, B, C, D, K] = lqrMatrices(theta0, psi0, phi0, u0, v0, w0, p0, q0, r0, T0, xi0, zeta0, M, g, rho2, R, H)
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

    
    % Calculate rotational inertias in body frame.
    Ixx = M*R^2/2;
    Iyy = M*(3*R^2+H^2)/12;
    Izz = Iyy;
    
    % Initialize state space matricies. Since they are sparse simply use zeros
    % and then input non-zero terms individually.
    % x = [x y z u v w phi theta psi p q r]'
    A = [
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 -q0 -p0 0 -g*cos(phi0)*sin(psi0)+g*cos(psi0)*sin(phi0)*sin(theta0) -g*cos(phi0)*cos(psi0)*cos(theta0) -g*cos(psi0)*sin(phi0)+g*cos(phi0)*sin(psi0)*sin(theta0) -v0 u0 0;
        0 0 0 0 r0 -q0 0 g*cos(psi0)*sin(theta0) g*cos(theta0)*sin(psi0) 0 -w0 v0;
        0 0 0 r0 0 p0 -g*sin(phi0)*sin(psi0)-g*cos(phi0)*cos(psi0)*sin(theta0) g*cos(psi0)*cos(theta0)*sin(phi0) g*cos(phi0)*cos(psi0)+g*sin(phi0)*sin(psi0)*sin(theta0) w0 0 -u0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 Iyy/Ixx*r0-Izz*r0/Ixx Iyy/Ixx*q0-Izz*q0/Ixx;
        0 0 0 0 0 0 0 0 0 -Ixx/Iyy*r0+Izz/Iyy*r0 0 -Ixx/Iyy*p0+Izz/Iyy*p0;
        0 0 0 0 0 0 0 0 0 Ixx*q0/Izz-Iyy*q0/Izz Ixx*p0/Izz-Iyy*p0/Izz 0;
    ];

    % u = [T tauR xi zeta]'
    B = [
        0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        cos(xi0)*cos(zeta0)/M 0 -T0*cos(zeta0)*sin(xi0)/M -T0*cos(xi0)*sin(zeta0)/M;
        sin(zeta0)/M 0 0 T0*cos(zeta0)/M;
        -cos(zeta0)*sin(xi0)/M 0 0 T0*sin(xi0)*sin(zeta0)/M -T0*cos(xi0)*cos(zeta0)/M;
        0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 1/Ixx 0 0;
        -rho2*cos(zeta0)*sin(xi0)/Iyy 0 -T0*rho2*cos(xi0)*cos(zeta0)/Iyy T0*rho2*sin(xi0)*sin(zeta0)/Iyy;
        -rho2*sin(zeta0)/Izz 0 0 -T0*rho2*cos(zeta0)/Izz;
    ];

    f = [
        0;
        0;
        0;
        q0*w0-r0*v0-g*cos(psi0)*cos(theta0)+(T0*xi0*cos(zeta0)*sin(x0)+T0*zeta0*cos(xi0)*sin(zeta0))/M-g*psi0*cos(theta0)*sin(psi0)-g*theta0*cos(psi0)*sin(theta0);
        T0*zeta0*cos(zeta0)/M-p0*w0+r0*u0+g*cos(phi0)*sin(psi0)+g*psi0*cos(phi0)*cos(psi0)+g*phi0*sin(phi0)*sin(psi0)-g*cos(psi0)*sin(phi0)*sin(theta0)+g*phi0*cos(phi0)*cos(psi0)*sin(theta0)+g*theta0*cos(psi0)*cos(theta0)*sin(phi0)-g*psi0*sin(phi0)*sin(psi0)*sin(theta0);
        cos(zeta0)*sin(xi0)/M-p0*v0+q0*u0-T0*zeta0*sin(xi0)*sin(zeta0)/M-g*sin(phi0)*sin(psi0)+T0*xi0*cos(xi0)*cos(zeta0)/M+g*phi0*cos(phi0)*sin(psi0)+g*psi0*cos(psi0)*sin(phi0)-g*cos(phi0)*cos(psi0)*sin(theta0)+g*theta0*cos(phi0)*cos(psi0)*cos(theta0)-g*phi0*cos(psi0)*sin(phi0)*sin(theta0)-g*psi0*cos(phi0)*sin(psi0)*sin(theta0);
        0;
        0;
        0;
        -Iyy*q0*r0/Ixx+Izz*q0*r0/Ixx;
        Ixx/Iyy*p0*r0-Izz/Iyy*p0*r0+T0*rho2*xi0*cos(xi0)*cos(zeta0)/Iyy-T0*rho2*zeta0*sin(xi0)*sin(zeta0)/Iyy;
        -Ixx*p0*q0/Izz+Iyy*p0*q0/Izz+T0*rho2*zeta0*cos(zeta0)/Izz;
    ];
   
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
