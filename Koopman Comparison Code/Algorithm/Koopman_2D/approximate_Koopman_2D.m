function K = approximate_Koopman_2D(N, x_GL, Ds, dynamics, p)
% Finite dimensional approximation of the Koopman Operator and
% its eigen-decomposition for 2D systems. 
%
% Args:
%       N: number of Gauss-Lobatto points
%       x_GL: Gauss-Lobatto points, each column corresponds to a component
%       Ds: Differentiation matrices, a struct
%       p: inherent parameters of the model, a struct
%       dynamics: dynamics of the model, function handle
% Returns:
%       K: koopman approximation of vectorized version
%       V: eigenvector matrix of vectorized version
%       L: eigenvalue matrix of vectorized version

x1_GL = x_GL(:, 1);
x2_GL = x_GL(:, 2);

D1 = Ds.D1;
D2 = Ds.D2;

[X, Y] = ndgrid(x1_GL, x2_GL);
% [X, Y] = meshgrid(x1_GL, x2_GL);
% X = X';
% Y = Y';

[F1, F2] = dynamics(X, Y, p);
I = eye(N);
K = diag(F1(:)) * kron(I, D1) + diag(F2(:)) * kron(D2, I);

