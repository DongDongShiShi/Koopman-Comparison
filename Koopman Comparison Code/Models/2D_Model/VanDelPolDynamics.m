function [F1, F2] = VandelPolDynamics(X, Y, p)

F1 = Y;
F2 = -X  + p * (1-X.^2) * Y;

