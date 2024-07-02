function [F1, F2, F3] = KraichnanOrszagDynamics(X, Y, Z, p)

F1 = p.k1 * Y .* Z;
F2 = p.k2 * X .* Z;
F3 = p.k3 * (-2) * X .* Y;
