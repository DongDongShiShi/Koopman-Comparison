function [K] = approximate_Koopman_1D(x_GL, D, dynamics,p)
% Finite dimensional approximation of the Koopman Operator and
% 
%
% Args:
%       x_GL: Gauss-Lobatto points
%       D: Differentiation matrix
%       dynamics: dynamics of the model, function handle
%       p: inherent parameters of the model, a struct
% Returns:
%       K: koopman approximation


f = dynamics(x_GL,p);
K = diag(f) * D;
