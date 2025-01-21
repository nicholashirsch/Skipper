
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
%             theta, phi, psi: I -> U Tait-Bryan angles (pitch, roll, yaw.)
%    thetaDot, phiDot, psiDot: Angular velocities of Tait-Bryan angles.
% thetaDDot, phiDDot, psiDDot: Angular accelerations of Tait-Bryan angles.
%                           T: Thrust.
%                          xi: First gimbal angle.
%                        zeta: Second gimbal angle.
%               Ixx, Iyy, Izz: Rigid body rotational inertias for a cylinder.
%                        rho2: Displacement of Skipper's CT relative to its CM.
%                           M: Skipper's mass.
%                           g: Acceleration due to gravity.
%                        tauR: Reaction wheel torque magnitude.
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
%                            rho: Position.
%                              v: Velocity.
%                              a: Acceleration.
%             alpha, beta, gamma: Euler angles (extrinsic).
%    alphaDot, betaDot, gammaDot: Angular velocity components.
% alphaDDot, betaDDot, gammaDDot: Angular acceleration components.
%                          omega: Angular velocity.
%                       omegaDot: Angular accelerations.
%                              S: Euler angle to angular veloctiy transformation matrix.
%                              I: Rotational inertia (rigid cylinder about its axis of symmetry).
%                              H: Angular momentum.
%                             Fg: Force due gravity.
%                             Ft: Force due to thrust.
%                          taur : Torque due to reaction wheels. 

% STEP 1A: KINEMATICS - LINEAR.
rho  = [x; y; z];             % Frame: N/A.
Rho2 = [0; 0; -rho2];         % Frame: N/A.
v    = [xDot; yDot; zDot];    % Frame: I.
a    = [xDDot; yDDot; zDDot]; % Frame: I.

% STEP 2: KINEMATICS - ANGULAR.
alphaDot = phiDot*cos(theta)*cos(psi) - thetaDot*sin(psi); % Frame: I.
betaDot = thetaDot*cos(psi) + phiDot*cos(theta)*sin(psi);  % Frame: I.
gammaDot = psiDot - phiDot*sin(theta);                     % Frame: I.

alphaDDot = phiDDot*cos(theta)*cos(phi) + ...
    phiDot*(-thetaDot*sin(theta)*cos(psi)-psiDot*cos(theta)*sin(psi)) ...
    - thetaDDot*sin(psi) - thetaDot*psiDot*cos(psi);
betaDDot = thetaDDot*cos(psi) - thetaDot*psiDot*sin(psi) ...
    + phiDDot*cos(theta)*sin(psi) ...
    + phiDot*(-thetaDot*sin(theta)*sin(psi)+psiDot*cos(theta)*cos(psi));
gammaDDot = psiDDot - phiDDot*sin(theta) - phiDot*thetaDot*cos(theta);

omega = [alphaDot; betaDot; gammaDot]; % Frame: I.
omegaDot = [alphaDDot; betaDDot; gammaDDot];

I     = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];            % Frame: U.
H     = I*omega;                           % Frame: I.
HDot  = I*omegaDot; % Frame: I.

% STEP 3: KINETICS - FORCES AND TORQUES.
Fg    = [0; 0; -M*g]; % Frame: I.
Ft    = [0; 0; T];   % Frame: T.
taur  = [0; 0; tauR]; % Frame: U

% STEP 4: EULER'S FIRST LAW.
% Note: Yields ux, uy, uz translational acceleration equations in V.
firstLawLHS = Tu2i*Tt2u*Ft + Fg;
firstLawRHS = M*a;
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
secondLawLHS = Tu2i*Tt2u*cross(Rho2, Ft) + Tu2i*taur;
secondLawRHS = HDot;
secondLaw    = secondLawLHS == secondLawRHS;

alphaEqnLHS = secondLawLHS(1);
alphaEqnRHS  = secondLawRHS(1);
alphaEqn    = alphaEqnLHS == alphaEqnRHS;

betaEqnLHS = secondLawLHS(2);
betaEqnRHS = secondLawRHS(2);
betaEqn    = betaEqnLHS == betaEqnRHS;

gammaEqnLHS = secondLawLHS(3);
gammaEqnRHS = secondLawRHS(3);
gammaEqn    = gammaEqnLHS == gammaEqnRHS;


% Solve equations of motion.
nonLinSol = solve([xEqn, yEqn, zEqn, alphaEqn, betaEqn, gammaEqn], ...
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
alphaTaylorLHS = taylor(alphaEqnLHS, eomVars, eomBasePoints, Order=O);
alphaTaylorRHS = taylor(alphaEqnRHS, eomVars, eomBasePoints, Order=O);
betaTaylorLHS = taylor(betaEqnLHS, eomVars, eomBasePoints, Order=O);
betaTaylorRHS = taylor(betaEqnRHS, eomVars, eomBasePoints, Order=O);
gammaTaylorLHS = taylor(gammaEqnLHS, eomVars, eomBasePoints, Order=O);
gammaTaylorRHS = taylor(gammaEqnRHS, eomVars, eomBasePoints, Order=O);


% Solve linearized equations of motion.
xTaylor = xTaylorLHS == xTaylorRHS;
yTaylor = yTaylorRHS == yTaylorLHS;
zTaylor = zTaylorLHS == zTaylorRHS;
alphaTaylor = alphaTaylorLHS == alphaTaylorRHS;
betaTaylor = betaTaylorLHS == betaTaylorRHS;
gammaTaylor = gammaTaylorLHS == gammaTaylorRHS;

linSol = solve([xTaylor, yTaylor, zTaylor, alphaTaylor, betaTaylor, gammaTaylor], ...  
    [xDDot, yDDot, zDDot,  phiDDot, thetaDDot, psiDDot]) % Linear solution.