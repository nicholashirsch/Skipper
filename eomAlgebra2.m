
% ________________________________________________________________________
%
%                      SKIPPER EQUATIONS OF MOTION v2
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  N. HIRSCH
% DATE:     1/18/2024
%
% OVERVIEW: THESE ARE THE EQUATIONS OF MOTION FOR SKIPPER. THE ARE DEFINED
%           IN THE MATLAB SYMBOLIC ENVIRONMENT UNSOLVED FOR USE IN FURTHER
%           MATHEMATICAL MANIPULATION.

clear; clc;

% I: Inertial frame.
% U: Vehicle-fixed frame.
% T: Thruster-fixed frame.

% Set up the MATLAB symbolic environment for the EOM derived from Euler's
% 1st and 2nd laws for rigid bodies. The following variables are used:
syms x y z ...                     % Displacements.
        xDot yDot zDot ...         % Velocities in I.
        u v w ...                  % Velocities in U.
        uDot vDot wDot ...         % Accelerations.
        theta phi psi ...          % Pitch, roll, and yaw.
        thetaDot phiDot psiDot ... % Euler angle rates.
        p q r ...                  % Angular rates in U.
        pDot qDot rDot ...         % Angular rate derivatives.
        xi zeta ...                % Gimbal angles.
        M Ixx Iyy Izz ...          % Linear and rotational inertias.
        rho1 rho2 ...     
        rho3x rho3y rho3z ...      % Lever arm distances.
        g ...                      % Acceleration due to gravity.
        T tauR ...                 % Control forces.
        alpha beta ...             % Wind angles.
        L D Q ...                  % Aerodynamic forces.
        Px Py Pz...                % Peturbation forces
        real % Avoid 'conj' popping up when performing matrix operations, 
             % always place after the last uncommented variable in 'syms' 
             % initialization.

% Generate transformation matrices.
Ti2u = euler2rMatrix(phi, 1)*euler2rMatrix(theta, 2)*...
        euler2rMatrix(psi, 3); % I -> U.
Tu2i = Ti2u.'; % U -> I.
Tu2t = euler2rMatrix(zeta, 3)*euler2rMatrix(xi, 2); % U -> T.
Tt2u = Tu2t.'; % T -> U.
Tu2w = euler2rMatrix(alpha, 3)*euler2rMatrix(beta, 2); % U -> W.
Tw2u = Tu2w'; % W -> U.

% Kinematics.
omega = [p q r]'; % Angular velocity.
alpha = [pDot qDot rDot]'; % Angular acceleration.
 
nu = [u v w]'; % Velocity.
a = [uDot vDot wDot]' + cross(omega, nu); % Acceleration. Need TT because 
                                         % components of velocity are in U.

I = diag([Ixx Iyy Izz]); % Rotational inertia tensor.

H = I*omega; % Angular momentum. THIS ONLY WORKS 
             % BECAUSE I IS DIAGONAL!
HDot = I*alpha + cross(omega, H); % Angular momentum derivative, need TT. 
                                  % Same rule regarding matrix
                                  % instead of tensor product as with H.

% Kinetics. Aerodynamic and peturbation forces are disabled for now.
F_aero = [0 0 0]'; % [-D Q -L]'; % Defined in W.
F_grav = [-M*g 0 0]'; % Defined in I.
F_thrust = [T 0 0]'; % Defined in T.
F_peturb = [0 0 0]'; % [Px Py Pz]'; % Defined in U.

tau_RCS = [tauR 0 0]';

% Convert lever arm distances to positions in U.
d1 = [rho1 0 0]'; % Defined in U.
d2 = [-rho2 0 0]'; % Defined in U.
d3 = [rho3x rho3y rho3z]'; % Defined in U.

% EOM derived via Newton-Euler laws. Note that all quantities are converted
% to U using rotation matricies. To make the use of 'taylor' possible later
% both sides of each law are first defined seperately and then combined.

first_law_LHS = Tw2u*F_aero + Ti2u*F_grav + Tt2u*F_thrust + F_peturb;
first_law_RHS = M*a;
first_law = first_law_LHS == first_law_RHS;

second_law_LHS = cross(d1, Tw2u*F_aero) + cross(d2, Tt2u*F_thrust) + ...
        cross(d3, F_peturb) + tau_RCS;
second_law_RHS = HDot; 
second_law = second_law_LHS == second_law_RHS; % Where Q = CM.

% First order differential equations relating linear and angular
% displacement to the body-frame velocities. Needed for state space.
S = [
    1 sin(phi)*tan(theta) cos(phi)*tan(theta);
    0 cos(phi) -sin(phi);
    0 sin(phi)*sec(theta) cos(phi)*sec(theta)
    ]; % Angular rates to euler angle rates conversion.

avel_rel_LHS = [phiDot thetaDot psiDot]'; % Euler angle rates.
avel_rel_RHS = S*omega;
avel_rel = avel_rel_LHS == avel_rel_RHS; % Relation betwen euler angle 
                                         % and body-frame rates.

vel_rel_LHS = [xDot yDot zDot]'; % Body-frame velocity components.
vel_rel_RHS = nu - [r*y + q*z -p*z+r*x -q*x+p*y]';
vel_rel = vel_rel_LHS == vel_rel_RHS; % Relation betwen inertial  
                                      % and body-frame velocities.


% Solve the EOM using 'solve'.
nonlinear_sol = solve([first_law, second_law], ...
    [uDot vDot wDot pDot qDot rDot])

% Linearize the nonlinear EOM using 'taylor' and generic initial
% conditions. The second order terms have already been solved for and thus
% need not have base points. Same goes for constants. Aerodynamic and
% peturbation terms are neglected for now.
syms theta0 phi0 psi0 ...
        v0 u0 w0 ...
        p0 q0 r0 ...                      
        xi0 zeta0 ...                
        T0 tauR0 ... 
        x0 y0 z0 ...
        real % Avoid 'conj' popping up when performing matrix operations, 
             % always place after the last uncommented variable in 'syms' 
             % initialization.

order = 2; % Order of Taylor series. O > 2 is nonlinear!

eom_variables = [x y z theta psi phi u v w p q r T xi zeta tauR];
eom_base_points = [x0 y0 z0 theta0 psi0 phi0 u0 v0 w0 p0 q0 r0 ...
        T0 xi0 zeta0 tauR0];

first_law_LHS_linear = taylor(first_law_LHS, eom_variables, ...
        eom_base_points, Order=order);
first_law_RHS_linear = taylor(first_law_RHS, eom_variables, ...
        eom_base_points, Order=order);
second_law_LHS_linear = taylor(second_law_LHS, eom_variables, ...
        eom_base_points, Order=order);
second_law_RHS_linear = taylor(second_law_RHS, eom_variables, ...
        eom_base_points, Order=order);

avel_rel_RHS_linear = taylor(avel_rel_RHS, eom_variables, ...
        eom_base_points, Order=order);
vel_rel_RHS_linear = taylor(vel_rel_RHS, eom_variables, ...
        eom_base_points, Order=order);

first_law_linear = first_law_LHS_linear == first_law_RHS_linear;
second_law_linear = second_law_LHS_linear == second_law_RHS_linear;

avel_rel_linear = avel_rel_LHS == avel_rel_RHS_linear
vel_rel_linear = vel_rel_LHS == vel_rel_RHS_linear


% Solve the linearized EOM.
linear_sol = solve([first_law_linear, second_law_linear], ...
    [uDot vDot wDot pDot qDot rDot])
