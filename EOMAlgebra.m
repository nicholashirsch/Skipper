
% ________________________________________________________________________
%
%                        SKIPPER EQUATIONS OF MOTION
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  N. HIRSCH, L. SIEMENAS
% DATE:     12/03/2024
%
% OVERVIEW: THESE ARE THE EQUATIONS OF MOTION FOR SKIPPER. THE ARE DEFINED
%           IN THE MATLAB SYMBOLIC ENVIRONMENT UNSOLVED FOR USE IN FURTHER
%           MATHEMATICAL MANIPULATION.

clear; clc;


% Set up the MATLAB symbolic environment for the EOM derived from Euler's
% 1st and 2nd laws for rigid bodies. The following variables are used:
%                          I:  Inertial frame.
%                          V:  Vehicle-fixed frame.
%                         CM: Center of mass.
%                         CT: Center of thrust.
%
%                     x, y, z: Positions.
%            xDot, yDot, zDot: Velocities.
%         xDDot, yDDot, zDDot: Accelerations.
%             theta, phi, psi: Euler angles to go from I to V.
%    thetaDot, phiDot, psiDot: Angular velocities of Euler angles in I.
% thetaDDot, phiDDot, psiDDot: Angular accelerations of Euler angles in I.
%                     p, q, r: Angular velocities of vehicle in V.
%            pDot, qDot, rDot: Angular accelerations of vehicle in V.
%                           T: Thrust.
%                          xi: First gimbal axis.
%                        zeta: Second gimbal axis.
%               Ixx, Iyy, Izz: Rigid body rotational inertias for a cylinder.
%                        rho2: Displacement of Skipper's CT relative to its CM.
%                           M: Skipper's mass.
%                           g: Acceleration due to gravity.
%                        tauR: Reaction wheel torque.
syms x y z xDot yDot zDot xDDot yDDot zDDot ...
    theta thetaDot thetaDDot psi psiDot psiDDot phi phiDot phiDDot...
    T xi zeta M g Ixx Iyy Izz rho2 tauR ...
    x0 y0 z0 xDot0 yDot0 zDot0 xDDot0 yDDot0 zDDot0 ...
    theta0 thetaDot0 thetaDDot0 psi0 psiDot0 psiDDot0 phi0 phiDot0 phiDDot0...
    T0 xi0 zeta0 tauR0


% p, q, r, pDot, qDot, and rDot are intermediate variables which may be
% defined in terms of the Euler angles and their derivatives.
p = phiDot - psiDot*sin(theta);
q = thetaDot*cos(phi) + psiDot*sin(phi)*cos(theta);
r = psiDot*cos(phi)*cos(theta) - thetaDot*sin(phi);
pDot = phiDDot - psiDDot*sin(theta) - psiDot*cos(theta)*thetaDot;
qDot = -thetaDot*sin(phi)*phiDot + thetaDDot*cos(phi) + psiDDot*sin(phi)*cos(theta) ...
    + psiDot*cos(phi)*phiDot*cos(theta) - psiDot*sin(phi)*sin(theta)*thetaDot;
rDot = psiDDot*cos(phi)*cos(theta) - psiDot*sin(phi)*phiDot*cos(theta) ...
    - psiDot*cos(phi)*sin(theta)*thetaDot - thetaDDot*sin(phi) - thetaDot*cos(phi)*phiDot;


% Euler's 1st law; in x, y, z axis order.
eqn1  = T*cos(xi)*cos(zeta) + M*g*sin(theta) ...
    == xDDot*cos(theta)*cos(psi) + yDDot*cos(theta)*sin(psi) - zDDot*sin(theta);

% Splitting each equation into right hand side and left hand side to do
% taylor linear approximations.
eqn1LHS = T*cos(xi)*cos(zeta) + M*g*sin(theta);
eqn1RHS = xDDot*cos(theta)*cos(psi) + yDDot*cos(theta)*sin(psi) - zDDot*sin(theta);


eqn2  = -T*sin(xi) - M*g*sin(phi)*sin(theta) ...
    == xDDot*(cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)) ...
    + yDDot*(cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)) + zDDot*sin(phi)*cos(theta);

eqn2LHS = -T*sin(xi) - M*g*sin(phi)*sin(theta);
eqn2RHS = xDDot*(cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)) ...
    + yDDot*(cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)) + zDDot*sin(phi)*cos(theta);


eqn3  = T*cos(xi)*sin(zeta) - M*g*cos(phi)*cos(theta) ...
    == xDDot*(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)) ...
    + yDDot*(-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)) + zDDot*cos(phi)*cos(theta);

eqn3LHS = T*cos(xi)*sin(zeta) - M*g*cos(phi)*cos(theta);
eqn3RHS = xDDot*(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)) ...
    + yDDot*(-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)) + zDDot*cos(phi)*cos(theta);


% Euler's 2nd law; in p, q, r order rotation order.
eqn4  = tauR == Ixx*pDot + p*q*r*(Izz-Iyy);

eqn4LHS = tauR;  % pretty stooopid o_O
eqn4RHS = Ixx*pDot + p*q*r*(Izz-Iyy);

eqn5  = rho2*T*cos(xi)*sin(zeta) == Iyy*qDot + p*q*r*(Izz-Ixx);

eqn5LHS = rho2*T*cos(xi)*sin(zeta);
eqn5RHS = Iyy*qDot + p*q*r*(Izz-Ixx);

eqn6  = rho2*T*sin(xi) == Izz*rDot + p*q*r*(Iyy-Ixx);

eqn6LHS = rho2*T*sin(xi);
eqn6RHS = Izz*rDot + p*q*r*(Iyy-Ixx);

% Solve equations of motion.
solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], ...
    [xDDot, yDDot, zDDot, thetaDDot psiDDot phiDDot])


% Spits out all the taylor series for all right and left hand sides.
taylor(eqn1LHS, [T, xi, zeta, theta], [T0, xi0, zeta0, theta0], Order=2)
taylor(eqn1RHS, [xDDot, theta, psi, yDDot], ...
    [xDDot0, theta0, psi0, yDDot0], Order=2)
taylor(eqn2LHS, [T, xi, phi, theta], [T0, xi0, phi0, theta0], Order=2)
taylor(eqn2RHS, [xDDot, yDDot, zDDot, phi, psi, theta], ...
    [xDDot0, yDDot0, zDDot0, phi0, psi0, theta0], Order=2)
taylor(eqn3LHS, [T, xi, zeta, phi, theta], ...
    [T0, xi0, zeta0, phi0, theta0], Order=2)
taylor(eqn3RHS, [xDDot, yDDot, zDDot, phi, psi, theta], ...
    [xDDot0, yDDot0, zDDot0, phi0, psi0, theta0], Order=2)
taylor(eqn4LHS, tauR, tauR0, Order=2)
taylor(eqn4RHS, [phiDDot, psiDDot, theta, phi, psi, thetaDot, ...
    phiDot, psiDot], [phiDDot0, psiDDot0, theta0, phi0, psi0, thetaDot0, ...
    phiDot0, psiDot0], Order=2)
taylor(eqn5LHS, [T, xi, zeta], [T0, xi0, zeta0], Order=2)
taylor(eqn5RHS, [thetaDDot, psiDDot, phi, psi, theta, phiDot, ...
    psiDot, thetaDot], [thetaDDot0, psiDDot0, phi0, psi0, theta0, ...
    phiDot0, psiDot0, thetaDot0], Order=2)
taylor(eqn6LHS, [T, xi], [T0, xi0], Order=2)
taylor(eqn6RHS, [thetaDDot, psiDDot, phi, psi, theta, phiDot, ...
    psiDot, thetaDot], [thetaDDot0, psiDDot0, phi0, psi0, theta0, ...
    phiDot0, psiDot0, thetaDot0], Order=2)