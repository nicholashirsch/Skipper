
% ________________________________________________________________________
%
%                          DYNAMICS SIMULATIONS
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  N. HIRSCH
% DATE:     12/04/2024
%
% OVERVIEW: EXECUTIVE FILE TO RUN SIMULATIONS USING EITHER THE NONLINEAR OR
%           LINEARIZED EQUATIONS OF MOTION FOR SKIPPER.

clear; clc; clf; close all;


% Define constants for the dynamic system.
%            CM: Center of mass.
%            CT: Center of thrust.
%
%             M: Mass.
%             g: Acceleration due to gravity.
%          rho2: Distance between the CM and CT.
% Ixx, Iyy, Izz: Rotational inertia.
M = 1; g = 10; rho2 = 1;
Ixx = 1; Iyy = 1; Izz = 1;

constants = [M g rho2 Izz Iyy Izz];


% Define equations of the controls as well as the base points to linearize
% the system about.
%   xi: First gimbal angle.
% zeta: Second gimbal angle.
%    T: Thrust.
% tauR: Reaction-wheel torque.
zeta = 0;
xi   = 0; 
T    = 15; 
tauR = 0;

T0 = 0.1; tauR0 = 0;

control    = [zeta xi T tauR];
basePoints = [T0 tauR0];
    

% Integrate linearized equations of motion using 'ode113'. p0 represents the
% starting point of the solver and is just the zero vector. 'ode113'
% outputs the discrete time grid it solved at as well as the phase, 'p' of
% the system which corresponds to its state.
p0 = [0 0 20 0 0 0 0 0 0 0 0 0]';
[t, p] = ode113(@(t, p)skipperODE(t, p, ...
    constants, control, basePoints), 0:0.01:10, p0);


% Extract Tait-Bryan states from the phase vector and use them to calculate
% extrinsic Euler angles - alpha, beta, and gamma - about the x, y, and z-
% axes respectively.
phi       = p(:, 7) ;
theta     = p(:, 8); 
psi       = p(:, 9) ;
phiDot    = p(:, 10);
thetaDot  = p(:, 11);
psiDot    = p(:, 12);
phiDDot   = (tauR0*p(8))/Ixx;
thetaDDot = -(tauR0*p(7))/Iyy;
psiDDot   = tauR/Izz;

alphaDot = phiDot.*cos(theta).*cos(psi) - thetaDot.*sin(psi);
betaDot  = thetaDot.*cos(psi) + phiDot.*cos(theta).*sin(psi);
gammaDot = psiDot - phiDot.*sin(theta);                      
alpha    = cumtrapz(alphaDot);
beta     = cumtrapz(betaDot);
gamma    = cumtrapz(gammaDot);


% Euler angles are on [0, 2pi] so apply the 'mod' function to restrict
% their domains.
alpha = mod(alpha, 2*pi);
beta  = mod(beta, 2*pi);
gamma = mod(gamma, 2*pi);


% Plot results.
subplot(2, 2, 1)
hold on

plot(t, p(:, 1))
plot(t, p(:, 2))
plot(t, p(:, 3))

xl = xlabel('Time, $t$');
yl = ylabel('Position, $[x, y, z]$');
ll = legend('$x(t)$', '$y(t)$', '$z(t)$', 'Location', 'southeast');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(ll,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;  
hold off

subplot(2, 2, 2)
hold on
plot(t, p(:, 4))
plot(t, p(:, 5))
plot(t, p(:, 6))

xl = xlabel('Time, $t$');
yl = ylabel('Velocity, $[\dot{x}, \dot{y}, \dot{z}]$');
ll = legend('$\dot{x}(t)$', '$\dot{y}(t)$', '$\dot{z}(t)$', 'Location', 'southeast');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(ll,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;       
hold off

subplot(2, 2, 3)
hold on
plot(t, alpha)
plot(t, beta)
plot(t, gamma)

xl = xlabel('Time, $t$');
yl = ylabel('Euler Angles, $[\alpha, \beta, \gamma]$');
ll = legend('$\alpha(t)$', '$\beta(t)$', '$\gamma(t)$', 'Location', 'southeast');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(ll,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;     
hold off

subplot(2, 2, 4)
hold on
plot(t, alphaDot)
plot(t, betaDot)
plot(t, gammaDot)

xl = xlabel('Time, $t$');
yl = ylabel('Euler Anglular Velocities, $[\dot{\alpha}, \dot{\beta}, \dot{\gamma}]$');
ll = legend('$\dot{\alpha}(t)$', '$\dot{\beta}(t)$', '$\dot{\gamma}(t)$', 'Location', 'southeast');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(ll,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;       
hold off


% The ODE function is defnied locally for convience.
function pDot = skipperODE(t, p, constants, control, basePoints)
% ________________________________________________________________________
%
%                                skipperODE.m
% ________________________________________________________________________
% 
% BY:       FLORIDA ROCKET LAB - GNC SUBTEAM
% EDITORS:  N. HIRSCH
% DATE:     12/07/2024
%
% OVERVIEW: ODE OBJECT WHICH CONTAINS THE EQUATIONS OF MOTION OF SKIPPER IN
%           FIRST ORDER FORM.
%
% INPUTS:
%              t: Discrete time points chosen by 'ode113' to solve system at.
%              p: Phase vector containing all states of the system
%      constants: Unchanging values of the system.
%        control: External inputs into the system.
%     basePoints: Base points of Taylor series used to linearize system.
% OUTPUTS:
%           pDot: Derivative of the phase. Internal variable passed to
%                 'ode113'.
    
    % Extract from inputs.
    M     = constants(1);
    g     = constants(2);
    % rho2  = constants(3); Not used.
    Ixx   = constants(4);
    Iyy   = constants(5);
    Izz   = constants(6);
     
    zeta  = control(1);
    xi    = control(2);
    T     = control(3);
    tauR  = control(4);
    
    T0    = basePoints(1);
    tauR0 = basePoints(2);
    

    % Initialize pDot and then define each of its components.
    pDot    = zeros(12,1);

    pDot(1)  = p(4);                   % xDot.
    pDot(2)  = p(5);                   % yDot.
    pDot(3)  = p(6);                   % zDot.
    pDot(4)  = (T0*p(8) + T0*xi)/M;    % xDDot.
    pDot(5)  = -(T0*p(7) + T0*zeta)/M; % yDDot.
    pDot(6)  = (T - M*g)/M;            % zDDot.
    pDot(7)  = p(10);                  % phiDot.
    pDot(8)  = p(11);                  % thetaDot.
    pDot(9)  = p(12);                  % psiDot.
    pDot(10) = (tauR0*p(8))/Ixx;       % phiDDot.
    pDot(11) = -(tauR0*p(7))/Iyy;      % thetaDDot.
    pDot(12) = tauR/Izz;               % psiDDot.
end