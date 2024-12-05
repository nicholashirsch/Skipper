
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
%                          I: Inertial frame.
%                          U: Vehicle-fixed frame.
%                          T: Thruster-fixed frame.
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
%    x0 y0 z0 xDot0 yDot0 zDot0 xDDot0 yDDot0 zDDot0 ...
%    theta0 thetaDot0 thetaDDot0 psi0 psiDot0 psiDDot0 phi0 phiDot0 phiDDot0...
%    T0 xi0 zeta0 tauR0


% Define rotation matricies used in shifting all terms in the EOM to V.
Tu2t = euler2rMatrix(zeta, 2)*euler2rMatrix(xi, 3); % Thrust to vehicle.
Tt2u = Tu2t.';

Ti2u = euler2rMatrix(phi, 1)*euler2rMatrix(theta, 2)*...
    euler2rMatrix(psi, 3); % Vehicle to inertial.
Tu2i = Ti2u.';


% The equations of motion of Skipper can be derived via fours steps:
%   1: Linear kinematics.
%   2: Angular kinematics.
%   3: Define forces and torques.
%   4: Euler's 1st law of rigid bodies.
%   5: Euler's 2nd law of rigid bodies
% In addition, the following variables are used:
%   rho: Position.
%     v: Velocity.
%     a: Acceleration.
% omega: Angular velocity.
%     I: Rotational inertia (rigid cylinder about its axis of symmetry).
%     H: Angular momentum.
%    Fg: Force due gravity.
%    Ft: Force due to thrust.
% taur : Torque due to reaction wheels. 

% STEP 1A: KINEMATICS - LINEAR.
rho  = [x; y; z];             % Frame: N/A.
rho2 = [-rho2; 0; 0];         % Frame: N/A.
v    = [xDot; yDot; zDot];    % Frame: I.
a    = [xDDot; yDDot; zDDot]; % Frame: I.

% STEP 2: KINEMATICS - ANGULAR.
% Note: p, q, r, pDot, qDot, and rDot are intermediate variables which may be
%       defined in terms of the Euler angles and their derivatives.
p    = phiDot - psiDot*sin(theta);
q    = thetaDot*cos(phi) + psiDot*sin(phi)*cos(theta);
r    = psiDot*cos(phi)*cos(theta) - thetaDot*sin(phi);
pDot = phiDDot - psiDDot*sin(theta) - psiDot*cos(theta)*thetaDot;
qDot = -thetaDot*sin(phi)*phiDot + thetaDDot*cos(phi) + psiDDot*sin(phi)*cos(theta) ...
        + psiDot*cos(phi)*phiDot*cos(theta) - psiDot*sin(phi)*sin(theta)*thetaDot;
rDot = psiDDot*cos(phi)*cos(theta) - psiDot*sin(phi)*phiDot*cos(theta) ...
        - psiDot*cos(phi)*sin(theta)*thetaDot - thetaDDot*sin(phi) ...
        - thetaDot*cos(phi)*phiDot;

omega = [p; q; r];                    % Frames: I -> U.
I     = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];  % Frame:  V.
H     = I*omega;                      % Frame:  I.
HDot  = [Ixx*pDot + p*q*r*(Izz-Iyy);  
         Iyy*qDot + p*q*r*(Izz-Ixx);  % Frame:  I.
         Izz*rDot + p*q*r*(Iyy-Ixx)];

% STEP 3: KINETICS - FORCES AND TORQUES.
Fg    = [0; 0; -M*g]; % Frame: I.
Ft    = [T; 0; 0];    % Frame: T.
taur  = [tauR; 0; 0]; % Frame: U

% STEP 4: EULER'S FIRST LAW.
% Note: Yields ux, uy, uz translational acceleration equations in V.
firstLawLHS = Tt2u*Ft + Ti2u*Fg;
firstLawRHS = Ti2u*M*a;

xEqnLHS = firstLawLHS(1);
xEqnRHS = firstLawRHS(1);
xEqn    = xEqnLHS == xEqnRHS;

yEqnLHS = firstLawLHS(2);
yEqnRHS = firstLawRHS(2);
yEqn    = yEqnLHS == yEqnRHS;

zEqnLHS = firstLawLHS(3);
zEqnRHS = firstLawRHS(3);
zEqn    = zEqnLHS == zEqnRHS;

% STEP 5: EULER'S SECOND LAW.
% Note: Yields ux, uy, and uz rotational acceleration equations in V.
secondLawLHS = Tt2u*cross(rho2, Ft) + taur;
secondLawRHS = HDot;

pEqnLHS = secondLawLHS(1);
pEqnRHS = secondLawRHS(1);
pEqn    = pEqnLHS == pEqnRHS;

qEqnLHS = secondLawLHS(2);
qEqnRHS = secondLawRHS(2);
qEqn    = qEqnLHS == qEqnRHS;

rEqnLHS = secondLawLHS(3);
rEqnRHS = secondLawRHS(3);
rEqn    = rEqnLHS == rEqnRHS;

% Euler's 1st law; in x, y, z axis order.
%eqn1LHS = T*cos(xi)*cos(zeta) + M*g*sin(theta);
%eqn1RHS = xDDot*cos(theta)*cos(psi) + yDDot*cos(theta)*sin(psi) - zDDot*sin(theta);
%eqn1    = eqn1LHS == eqn1RHS;
%
%eqn2LHS = -T*sin(xi) - M*g*sin(phi)*sin(theta);
%eqn2RHS = xDDot*(cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)) ...
%    + yDDot*(cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)) + zDDot*sin(phi)*cos(theta);
%eqn2    = eqn2LHS == eqn2RHS;
%
%eqn3LHS = T*cos(xi)*sin(zeta) - M*g*cos(phi)*cos(theta);
%eqn3RHS = xDDot*(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)) ...
%    + yDDot*(-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)) + zDDot*cos(phi)*cos(theta);
%eqn3    = eqn3LHS == eqn3RHS;
%
%
% Euler's 2nd law; in p, q, r order rotation order.
%eqn4LHS = tauR;  % pretty stooopid o_O
%eqn4RHS = Ixx*pDot + p*q*r*(Izz-Iyy);
%eqn4    = eqn4LHS == eqn4RHS;
%
%eqn5LHS = rho2*T*cos(xi)*sin(zeta);
%eqn5RHS = Iyy*qDot + p*q*r*(Izz-Ixx);
%eqn5    = eqn5LHS == eqn5RHS;
%
%eqn6LHS = rho2*T*sin(xi);
%eqn6RHS = Izz*rDot + p*q*r*(Iyy-Ixx);
%eqn6    = eqn6LHS == eqn6RHS;


% Solve equations of motion.
nonLinSol = solve([xEqn, yEqn, zEqn, pEqn, qEqn, rEqn], ...
    [xDDot, yDDot, zDDot, thetaDDot psiDDot phiDDot]) % Nonlinear solution.


% Spits out all the taylor series for all right and left hand sides. Define
% base points for approximations here. Some rules:
%   - By default all angles should be approximated about 0 [deg], this is
%     accurate till roughly 30 [deg].
%   - Approximation fails if T0 = 0, reason for this is unknown.
O = 2; % Order of Taylor series. Note that O > 2 is nonlinear!

    x0 = 0;      y0 = 0;      z0 = 0;
 xDot0 = 0;   yDot0 = 0;   zDot0 =  0;
xDDot0 = 0;  yDDot0 = 0;  zDDot0 = 0;

theta0 = 0;  thetaDot0 = 0;  thetaDDot0 = 0;
  psi0 = 0;    psiDot0 = 0;    psiDDot0 = 0;
  phi0 = 0;    phiDot0 = 0;    phiDDot0 = 0;
  xi0  = 0; 
 zeta0 = 0;

%syms T0 tauR0 % Useful for seeing where there show up in linearized solutions. 
              % Comment out for a specific solutions.
T0 = 1; tauR0 = 0;

xTaylorLHS = taylor(xEqnLHS, [T, xi, zeta, theta], [T0, xi0, zeta0, theta0], Order=O);
xTaylorRHS = taylor(xEqnRHS, [xDDot, theta, psi, yDDot], ...
    [xDDot0, theta0, psi0, yDDot0], Order=O);
yTaylorLHS = taylor(yEqnLHS, [T, xi, phi, theta], [T0, xi0, phi0, theta0], Order=O);
yTaylorRHS = taylor(yEqnRHS, [xDDot, yDDot, zDDot, phi, psi, theta], ...
    [xDDot0, yDDot0, zDDot0, phi0, psi0, theta0], Order=O);
zTaylorLHS = taylor(zEqnLHS, [T, xi, zeta, phi, theta], ...
    [T0, xi0, zeta0, phi0, theta0], Order=O);
zTaylorRHS = taylor(zEqnRHS, [xDDot, yDDot, zDDot, phi, psi, theta], ...
    [xDDot0, yDDot0, zDDot0, phi0, psi0, theta0], Order=O);
pTaylorLHS = taylor(pEqnLHS, tauR, tauR0, Order=O);
pTaylorRHS = taylor(pEqnRHS, [phiDDot, psiDDot, theta, phi, psi, thetaDot, ...
    phiDot, psiDot], [phiDDot0, psiDDot0, theta0, phi0, psi0, thetaDot0, ...
    phiDot0, psiDot0], Order=O);
qTaylorLHS = taylor(qEqnLHS, [T, xi, zeta], [T0, xi0, zeta0], Order=O);
qTaylorRHS = taylor(qEqnRHS, [thetaDDot, psiDDot, phi, psi, theta, phiDot, ...
    psiDot, thetaDot], [thetaDDot0, psiDDot0, phi0, psi0, theta0, ...
    phiDot0, psiDot0, thetaDot0], Order=O);
rTaylorLHS = taylor(rEqnLHS, [T, xi], [T0, xi0], Order=O);
rTaylorRHS = taylor(rEqnRHS, [thetaDDot, psiDDot, phi, psi, theta, phiDot, ...
    psiDot, thetaDot], [thetaDDot0, psiDDot0, phi0, psi0, theta0, ...
    phiDot0, psiDot0, thetaDot0], Order=O);


% Solve linearized equations of motion.
xTaylor = xTaylorLHS == xTaylorRHS;
yTaylor = yTaylorRHS == yTaylorLHS;
zTaylor = zTaylorLHS == zTaylorRHS;
pTaylor = pTaylorLHS == pTaylorRHS;
qTaylor = qTaylorLHS == qTaylorRHS;
rTaylor = rTaylorLHS == rTaylorRHS;

linSol = solve([xTaylor, yTaylor, zTaylor, pTaylor, qTaylor, rTaylor], ...
    [xDDot, yDDot, zDDot, thetaDDot psiDDot phiDDot]) % Linear solution.