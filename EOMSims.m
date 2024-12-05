
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


% Define linearized equations of motion.
T0    = 0;
tauR0 = 0;

p0 = [0 0 20 0 0 0 0 0 0 0 0 0]';

[t, p] = ode113(@(t, p)skipperODE(t, p), [0, 25], p0);


figure(1)
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

figure(2)
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

figure(3)
hold on
plot(t, p(:, 7))
plot(t, p(:, 8))
plot(t, p(:, 9))

xl = xlabel('Time, $t$');
yl = ylabel('Tait-Bryan Angles, $[\phi, \theta, \psi]$');
ll = legend('$\phi(t)$', '$\theta(t)$', '$\psi(t)$', 'Location', 'southeast');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(ll,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;       
hold off

figure(4)
hold on
plot(t, p(:, 10))
plot(t, p(:, 11))
plot(t, p(:, 12))

xl = xlabel('Time, $t$');
yl = ylabel('Tait-Bryan Anglular Velocities, $[\dot{\phi}, \dot{\theta}, \dot{\psi}]$');
ll = legend('$\dot{\phi}(t)$', '$\dot{\theta}(t)$', '$\dot{\psi}(t)$', 'Location', 'southeast');
set(xl,'FontSize',16,'Interpreter','LaTeX');
set(yl,'FontSize',16,'Interpreter','LaTeX');
set(ll,'FontSize',16,'Interpreter','LaTeX');
set(gca,'FontSize',16,'TickLabelInterpreter','LaTeX');
grid on;       
hold off

function pDot = skipperODE(t, p)
    
    M = 1;
    zeta = 0;
    xi = 0;
    T = 0;
    g = 10;
    Ixx = 1;
    Iyy = 1;
    Izz = 1;
    tauR = 0;
    T0 = 1;
    tauR0 = 1;
    rho2 = 1;
    pDot    = zeros(12,1);

    pDot(1)  = p(4); % xDot.
    pDot(2)  = p(5); % yDot.
    pDot(3)  = p(6); % zDot.
    pDot(4)  = (T0*xi + M*g*p(8))/M;  % xDDot.
    pDot(5)  = -(T0*zeta + M*g*p(7))/M; % yDDot.
    pDot(6)  = (T - M*g)/M; % zDDot.
    pDot(7)  = p(10); % phiDot.
    pDot(8)  = p(11); % thetaDot.
    pDot(9)  = p(12); % psiDot.
    pDot(10) = 0; % phiDDot.
    pDot(11) = T*rho2/Iyy; % thetaDDot.
    pDot(12) = (tauR + T0*rho2*zeta)/Izz; % psiDDot.
end