function [K,y0] = KSLinearization(dynamics,p,r,x0,N)

%parameter default
op=1;

% space
lb = x0 - r;
ub = x0 + r;

% compute Gauss-Lobatto points and differentiation matrices
[x_GL, D_r] = compute_diffMat_1D(N, lb, ub, op);
y0=x_GL;

% compute finite dimensional Koopman approximation matrix K
K = approximate_Koopman_1D(x_GL, D_r, dynamics, p);
 
end

