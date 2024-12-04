


clear; clc;

syms x y z xDot yDot zDot xDDot yDDot zDDot ...
    theta thetaDot thetaDDot psi psiDot psiDDot phi phiDot phiDDot...
    T zeta xi M g Ixx Iyy Izz rho2 tauR

p = phiDot - psiDot*sin(theta);
q = thetaDot*cos(phi) + psiDot*sin(phi)*cos(theta);
r = psiDot*cos(phi)*cos(theta) - thetaDot*sin(phi);
pDot = phiDDot - psiDDot*sin(theta) - psiDot*cos(theta)*thetaDot;
qDot = -thetaDot*sin(phi)*phiDot + thetaDDot*cos(phi) + psiDDot*sin(phi)*cos(theta) + psiDot*cos(phi)*phiDot*cos(theta) - psiDot*sin(phi)*sin(theta)*thetaDot;
rDot = psiDDot*cos(phi)*cos(theta) - psiDot*sin(phi)*phiDot*cos(theta) - psiDot*cos(phi)*sin(theta)*thetaDot - thetaDDot*sin(phi) - thetaDot*cos(phi)*phiDot;

eqn1  = T*cos(xi)*cos(zeta) + M*g*sin(theta) == xDDot*cos(theta)*cos(psi) + yDDot*cos(theta)*sin(psi) - zDDot*sin(theta);
eqn2  = -T*sin(xi) - M*g*sin(phi)*sin(theta) == xDDot*(cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)) + yDDot*(cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)) + zDDot*sin(phi)*cos(theta);
eqn3  = T*cos(xi)*sin(zeta) - M*g*cos(phi)*cos(theta) == xDDot*(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)) + yDDot*(-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)) + zDDot*cos(phi)*cos(theta);
eqn4  = tauR == Ixx*pDot + p*q*r*(Izz-Iyy);
eqn5  = rho2*T*cos(xi)*sin(zeta) == Iyy*qDot + p*q*r*(Izz-Ixx);
eqn6  = rho2*T*sin(xi) == Izz*rDot + p*q*r*(Iyy-Ixx);

solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [xDDot, yDDot, zDDot, thetaDDot psiDDot phiDDot])