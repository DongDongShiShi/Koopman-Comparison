function [F1, F2] = LimitCycleDynamics(X, Y, p)

F1 =  - X - Y + p.k * X ./ sqrt(X.^2 + Y.^2);
F2 =    X - Y + p.k * Y ./ sqrt(X.^2 + Y.^2);
