
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
%             theta, phi, psi: I -> U euler angles (pitch, roll, yaw.)
%    thetaDot, phiDot, psiDot: Angular velocities of Euler angles in I.
% thetaDDot, phiDDot, psiDDot: Angular accelerations of Euler angles in I.
%                           T: Thrust.
%                          xi: First gimbal angle.
%                        zeta: Second gimbal angle.
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




% Define rotation matricies used in shifting all terms in the EOM to U.
% The following transformations are used:
%   I -> U: 3-2-1 transformation through psi, theta, phi.
%   U -> T: 2-1 transformation through xi, zeta.
Ti2u = euler2rMatrix(phi, 1)*euler2rMatrix(theta, 2)*euler2rMatrix(psi, 3); % I -> U.
Tu2i = Ti2u.';

Tu2t = euler2rMatrix(zeta, 1)*euler2rMatrix(xi, 2); % U -> T.
Tt2u = Tu2t.';


% The equations of motion of Skipper can be derived via fours steps:
%   1: Linear kinematics.
%   2: Angular kinematics.
%   3: Define forces and torques.
%   4: Euler's 1st law of rigid bodies.
%   5: Euler's 2nd law of rigid bodies
% In addition, the following variables are used:
%              rho: Position.
%                v: Velocity.
%                a: Acceleration.
%          p, q, r: Angular velocity components.
% pDot, qDot, rDot: Angular acceleration components.
%            omega: Angular velocity.
%         omegaDot: Angular accelerations.
%                S: Euler angle to angular veloctiy transformation matrix.
%                I: Rotational inertia (rigid cylinder about its axis of symmetry).
%                H: Angular momentum.
%               Fg: Force due gravity.
%               Ft: Force due to thrust.
%            tauR : Torque due to reaction wheels. 

% STEP 1A: KINEMATICS - LINEAR.
rho  = [x; y; z];             % Frame: N/A.
rho2 = [-rho2; 0; 0];         % Frame: N/A.
v    = [xDot; yDot; zDot];    % Frame: I.
a    = [xDDot; yDDot; zDDot]; % Frame: I.

% STEP 2: KINEMATICS - ANGULAR.
p = phiDot - psiDot*sin(theta);                     % Frame: I.
q = thetaDot*cos(phi) + psiDot*sin(phi)*cos(theta); % Frame: I.
r = psiDot*cos(phi)*cos(theta) - thetaDot*sin(phi); % Frame: I.

pDot = phiDDot - psiDDot*sin(theta) - psiDot*thetaDot*cos(theta);
qDot = thetaDDot*cos(phi) - thetaDot*phiDot*sin(phi) +...
        psiDDot*sin(theta)*cos(theta) + ...
        psiDot*(phiDot*cos(phi)*cos(theta)-thetaDot*cos(phi)*sin(theta));
rDot = psiDDot*cos(phi)*cos(theta) + ...
    psiDot*(-phiDot*sin(phi)*cos(theta)-thetaDot*cos(phi)*sin(theta)) ...
    - thetaDDot*sin(phi) - thetaDot*phiDot*cos(phi);

omega = [p; q; r]; % Frame: I.
omegaDot = [pDot; qDot; rDot];

I     = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];  % Frame: U.
H     = I*omega;                      % Frame: I.
HDot  = I*omegaDot;                   % Frame: I.

% STEP 3: KINETICS - FORCES AND TORQUES.
Fg    = [0; 0; -M*g]; % Frame: I.
Ft    = [0; 0; T];   % Frame: T.
taur  = [0; 0; tauR]; % Frame: U

% STEP 4: EULER'S FIRST LAW.
% Note: Yields ux, uy, uz translational acceleration equations in V.
firstLawLHS = Tt2u*Ft + Ti2u*Fg;
firstLawRHS = Ti2u*M*a;
firstLaw    = firstLawLHS == firstLawRHS;

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
secondLawRHS = Ti2u*HDot;
secondLaw    = secondLawLHS == secondLawRHS;

pEqnLHS = secondLawLHS(1);
pEqnRHS = secondLawRHS(1);
pEqn    = pEqnLHS == pEqnRHS;

qEqnLHS = secondLawLHS(2);
qEqnRHS = secondLawRHS(2);
qEqn    = qEqnLHS == qEqnRHS;

rEqnLHS = secondLawLHS(3);
rEqnRHS = secondLawRHS(3);
rEqn    = rEqnLHS == rEqnRHS;


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
 xDot0 = 0;   yDot0 = 0;   zDot0 = 0;
xDDot0 = 0;  yDDot0 = 0;  zDDot0 = 0;

theta0 = 0;  thetaDot0 = 0;  thetaDDot0 = 0;
  psi0 = 0;    psiDot0 = 0;    psiDDot0 = 0;
  phi0 = 0;    phiDot0 = 0;    phiDDot0 = 0;
  xi0  = 0; 
 zeta0 = 0;

syms T0 tauR0 % Useful for seeing where there show up in linearized solutions. 
              % Comment out for a specific solutions.
%T0 = 1; tauR0 = 1;

eomVars = [x y z xDot yDot zDot xDDot yDDot zDDot theta thetaDot thetaDDot ...
    psi psiDot psiDDot phi phiDot phiDDot T xi zeta tauR];
eomBasePoints = [x0 y0 z0 xDot0 yDot0 zDot0 xDDot0 yDDot0 zDDot0 theta0 ...
    thetaDot0 thetaDDot0 psi0 psiDot0 psiDDot0 phi0 phiDot0 phiDDot0 T0 xi0 zeta0 tauR0];

xTaylorLHS = taylor(xEqnLHS, eomVars, eomBasePoints, Order=O);
xTaylorRHS = taylor(xEqnRHS, eomVars, eomBasePoints, Order=O);
yTaylorLHS = taylor(yEqnLHS, eomVars, eomBasePoints, Order=O);
yTaylorRHS = taylor(yEqnRHS, eomVars, eomBasePoints, Order=O);
zTaylorLHS = taylor(zEqnLHS, eomVars, eomBasePoints, Order=O);
zTaylorRHS = taylor(zEqnRHS, eomVars, eomBasePoints, Order=O);
pTaylorLHS = taylor(pEqnLHS, eomVars, eomBasePoints, Order=O);
pTaylorRHS = taylor(pEqnRHS, eomVars, eomBasePoints, Order=O);
qTaylorLHS = taylor(qEqnLHS, eomVars, eomBasePoints, Order=O);
qTaylorRHS = taylor(qEqnRHS, eomVars, eomBasePoints, Order=O);
rTaylorLHS = taylor(rEqnLHS, eomVars, eomBasePoints, Order=O);
rTaylorRHS = taylor(rEqnRHS, eomVars, eomBasePoints, Order=O);


% Solve linearized equations of motion.
xTaylor = xTaylorLHS == xTaylorRHS;
yTaylor = yTaylorRHS == yTaylorLHS;
zTaylor = zTaylorLHS == zTaylorRHS;
pTaylor = pTaylorLHS == pTaylorRHS;
qTaylor = qTaylorLHS == qTaylorRHS;
rTaylor = rTaylorLHS == rTaylorRHS;

linSol = solve([xTaylor, yTaylor, zTaylor, pTaylor, qTaylor, rTaylor], ...  
    [xDDot, yDDot, zDDot,  phiDDot, thetaDDot, psiDDot]) % Linear solution.