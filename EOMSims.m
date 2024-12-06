
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

    M = 1;
    zeta = 0;
    xi = 0;
    T = 15;
    g = 10;
    Ixx = 1;
    Iyy = 1;
    Izz = 1;
    tauR = 0;
    T0 = 0.1;
    tauR0 = 0;
    rho2 = 1;

p0 = [0 0 20 0 0 0 0 0 0 0 0 0]';

[t, p] = ode113(@(t, p)skipperODE(t, p, M, zeta, xi, T, g, Ixx, Iyy, Izz, tauR, T0, tauR0), 0:0.01:10, p0);

phi = p(:, 7) ;
theta = p(:, 8); 
psi = p(:, 9) ;
phiDot =p(:, 10);
thetaDot = p(:, 11);
psiDot = p(:, 12);
phiDDot = (tauR0*p(8))/Ixx; % phiDDot.
thetaDDot = -(tauR0*p(7))/Iyy; % thetaDDot.
psiDDot = tauR/Izz; % psiDDot.

alphaDot = phiDot.*cos(theta).*cos(psi) - thetaDot.*sin(psi); % Frame: I.
betaDot = thetaDot.*cos(psi) + phiDot.*cos(theta).*sin(psi);  % Frame: I.
gammaDot = psiDot - phiDot.*sin(theta);                       % Frame: I.

alpha = cumtrapz(alphaDot);
beta = cumtrapz(betaDot);
gamma = cumtrapz(gammaDot);

alpha = mod(alpha, 2*pi);
beta = mod(beta, 2*pi);
gamma = mod(gamma, 2*pi);

hold on
subplot(2, 2, 1)

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

function pDot = skipperODE(t, p, M, zeta, xi, T, g, Ixx, Iyy, Izz, tauR, T0, tauR0)
    pDot    = zeros(12,1);

    pDot(1)  = p(4); % xDot.
    pDot(2)  = p(5); % yDot.
    pDot(3)  = p(6); % zDot.
    pDot(4)  = (T0*p(8) + T0*xi)/M;  % xDDot.
    pDot(5)  = -(T0*p(7) + T0*zeta)/M; % yDDot.
    pDot(6)  = (T - M*g)/M; % zDDot.
    pDot(7)  = p(10); % phiDot.
    pDot(8)  = p(11); % thetaDot.
    pDot(9)  = p(12); % psiDot.
    pDot(10) = (tauR0*p(8))/Ixx; % phiDDot.
    pDot(11) = -(tauR0*p(7))/Iyy; % thetaDDot.
    pDot(12) = tauR/Izz; % psiDDot.
end