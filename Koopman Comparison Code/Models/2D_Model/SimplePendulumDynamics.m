function [F1, F2] = SimplePendulumDynamics(X, Y, p)

F1 = Y;
F2 = - (p.g / p.L) * sin(X);
