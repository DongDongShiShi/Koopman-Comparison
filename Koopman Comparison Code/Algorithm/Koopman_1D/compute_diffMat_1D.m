function [x_GL, D] = compute_diffMat_1D(N, lb, ub, option)
% This function computes the Gauss-Lobatto points and the differentiation
% matrix as in the collocation method. 
%
% Args:
%       N: number of Gauss-Lobatto points
%       lb: lower bound
%       ub: upper bound
%       option: 1. Chebyshev-Gauss-Lobatto
%               2. Legendre-Gauss-Lobatto
% Return:
%       x_GL: Gauss-Lobatto points
%       D: Differentiation matix

n = N - 1;

if option == 1
    [x_GL, D] = chebyshevDiff(n, lb, ub);
elseif option == 2
    [x_GL, D] = legendreDiff(n, lb, ub);
else
    error('Option Value Error: option can only be 1 (Chebyshev) or 2 (Legendre) ...');
end
