function [K] = approximate_Koopman_3D(N, x_GL, Ds, dynamics, p)
% Finite dimensional approximation of the Koopman Operator and
% its eigen-decomposition for 3D systems. 
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
x3_GL = x_GL(:, 3);

D1 = Ds.D1;
D2 = Ds.D2;
D3 = Ds.D3;

[X, Y, Z] = ndgrid(x1_GL, x2_GL, x3_GL);
% [X, Y, Z] = meshgrid(x1_GL, x2_GL, x3_GL);
% X = pagetranspose(X);
% Y = pagetranspose(Y);


[F1, F2, F3] = dynamics(X, Y, Z, p);
I = eye(N);
K = diag(F1(:)) * kron(kron(I, I), D1) + diag(F2(:)) * kron(kron(I, D2), I) + diag(F3(:)) * kron(kron(D3, I), I);


