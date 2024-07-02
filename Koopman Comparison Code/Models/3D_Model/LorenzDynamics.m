function [F1, F2, F3] = LorenzDynamics(X, Y, Z, p)

F1 = p.alpha * (Y - X);
F2 = X .* (p.beta - Z) - Y;
F3 = X .* Y - p.gamma * Z;
