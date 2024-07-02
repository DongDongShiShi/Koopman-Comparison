function [x1n, x2n, x3n] = numerical_solve_3D(ct, lt, V, L, C)
% Numerically solve the equations using the eigentensors 
% and eigenvalues given by the Koopman approx. 
%
% In particular, 
%   x = sum_i c_i * v_i(x_init) * exp(lambda* t),
% where v_i is eigenfunction value given by the eigentensor.
%
% Args:
%       ct: current time
%       lt: last time (of eigen-decomposition)
%       V: matrix of eigenvectors
%       L: matrix of eigenvalues
%       C: matrix of coefficients
% Returns:
%       x1n: numerical solution of the first component
%       x2n: numerical solution of the second component
%       x3n: numerical solution of the third component

[nrow, ~] = size(V);
idx = (nrow + 1) / 2;      % middle point

% solve a linear dynamics
lambda = diag(L);
v = V(idx, :);
c1 = C(:, 1);
c2 = C(:, 2);
c3 = C(:, 3);
delta_t = ct - lt;

x1n = real( v * (c1 .* exp(lambda .* delta_t)) );
x2n = real( v * (c2 .* exp(lambda .* delta_t)) );
x3n = real( v * (c3 .* exp(lambda .* delta_t)) );



