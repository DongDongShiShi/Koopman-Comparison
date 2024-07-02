function [xT] = KoopmanLinearization_1D(dynamics, p,T,N, x0, r, option)
% The adaptive spectral Koopman method numerically solves autonomous
% dynamical systems, using the properties of the Koopman operator. This
% method also leverages ideas from the spectral collocation methods.
% This implemetation is based on one-dimensional systems. 
% 
% Args:
%       dynamics: the dynamics f of the system, a function handle
%       p: parameter of the dynamics, a structure
%       g: observable of dynamics (default x itself)
%       T: terminal time
%       N: number of Gauss-Lobatto points 
%       n: number of output points
%       x0: initial state
%       r: radius
%       option: different type of differentiation matrix
%               1: Chebyshev (default) 2: Legendre 3.FiniteDifference
%
% Return:
%       gxnT: numerical solution of the system at T with respect to observable g
       

%% parameter parsing
assert(isa(dynamics, 'function_handle'), 'TypeError: ');
%assert(isstruct(p), 'TypeError: p must be a struct ...');
assert(T > 0, 'ValueError: T must be positive ...');
assert((N >= 3) && (mod(N, 2) == 1), 'ValueError: N must be at least 3 and odd ...');
assert(numel(x0) == numel(r), 'AttributeError: x0 and r must have the same length ...');
assert(numel(x0) == 1, 'AttributeError: x0 must be a scalar ...');
assert((option == 1) || (option == 2), 'ValueError: option can be only 1 or 2 or 3 ...');
op = 1;


%% Algorithm
lb = x0 - r;
ub = x0 + r;

% compute Gauss-Lobatto points and differentiation matrices
[x_GL_0, D] = compute_diffMat_1D(N, -1, 1, op);
[x_GL, D_r] = rescale_diffMat_1D(x_GL_0, D, lb, ub);

% compute finite dimensional Koopman approximation matrix K
K = approximate_Koopman_1D(x_GL, D_r, dynamics, p);

% compute projection coefficients g
% 1D default g(x)=x, which means y_0 = x_GL here
% Therefore we obtain dy/dt=Ky with y_0=x_GL as inital point
% solve linear homogenous ode by y=e^(A*t)*y_0
yT= expm(K*T)*x_GL;
xT=yT((N+1)/2);



