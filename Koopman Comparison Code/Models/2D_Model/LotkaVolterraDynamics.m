function [F1, F2] = LotkaVolterraDynamics(X, Y, p)

F1 = 1.1 * X - 0.4 * X .* Y;
F2 = 0.1 * X .* Y - 0.4 * Y;

