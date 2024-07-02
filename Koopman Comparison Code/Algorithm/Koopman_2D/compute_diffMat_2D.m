function [x_GL, Ds] = compute_diffMat_2D(N, lb, ub, option)
% This function computes the Gauss-Lobatto points and the differentiation
% matrix as in the collocation method. 
%
% Args:
%       N: number of Gauss-Lobatto points
%       lb: lower bounds
%       ub: upper bounds
%       option: 1. Chebyshev-Gauss-Lobatto
%               2. Legendre-Gauss-Lobatto
% Return:
%       x_GL: Gauss-Lobatto points, each column corresponds to a component
%       Ds: Differentiation matices, a struct

n = N - 1;

if option == 1
    [x1_GL, D1] = chebyshevDiff(n, lb(1), ub(1));
    [x2_GL, D2] = chebyshevDiff(n, lb(2), ub(2));
elseif option == 2
    [x1_GL, D1] = legendreDiff(n, lb(1), ub(1));
    [x2_GL, D2] = legendreDiff(n, lb(2), ub(2));
else
    error('Option Value Error: option can only be 1 (Chebyshev) or 2 (Legendre) ...');
end

x_GL = [x1_GL, x2_GL];
Ds.D1 = D1;
Ds.D2 = D2;