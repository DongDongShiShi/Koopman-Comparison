function [K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N, x0, r,option)
% The adaptive spectral Koopman method numerically solves autonomous
% dynamical systems, using the properties of the Koopman operator. This
% method also leverages ideas from the spectral collocation methods.
% This implemetation is based on two-dimensional systems. 
% 
% Args:
%       dynamics: the dynamics f of the system, a function handle
%       p: inherent parameters of the system
%       n: number of check points
%       T: terminal time
%       N: number of Gauss-Lobatto points 
%       x0: initial state
%       r: radius
%       option: different type of polynomials for differentiation
%               1: Chebyshev (default) 2: Legendre
%
% Return:
%       xnT: numerical solution of the system at T
%       n_decomp: number of eigen-decompositions

%% parameter parsing
assert(isa(dynamics, 'function_handle'), 'TypeError: ');
%assert(isstruct(p), 'TypeError: p must be a struct ...');
assert(T > 0, 'ValueError: T must be positive ...');
assert((N >= 3) && (mod(N, 2) == 1), 'ValueError: N must be at least 3 and odd ...');
assert(numel(x0) == numel(r), 'AttributeError: x0 and r must have the same length ...');
assert((option == 1) || (option == 2), 'ValueError: option can be only 1 or 2 ...');
op = option;
x0 = reshape(x0, 2, 1);
op = 1;


lb = [x0(1)-r(1), x0(2)-r(2)];
ub = [x0(1)+r(1), x0(2)+r(2)];
% compute Gauss-Lobatto points and differentiation matrices
[x_GL_0, Ds_0] = compute_diffMat_2D(N, [-1, -1], [1, 1], op);
[x_GL, Ds] = rescale_diffMat_2D(x_GL_0, Ds_0, lb, ub);

% compute finite dimensional Koopman approximation
K = approximate_Koopman_2D(N, x_GL, Ds, dynamics, p);

% compute Gauss-Lobatto points and differentiation matrices
[x_GL, D_r] = compute_diffMat_2D(N, lb, ub, op);

% compute finite dimensional Koopman approximation matrix K
K = approximate_Koopman_2D(N,x_GL, D_r, dynamics, p);
% compute projection coefficients g_0
% Now default g(x)=[x1;x2]
[init1,init2]=compute_initial_2D(x_GL);


