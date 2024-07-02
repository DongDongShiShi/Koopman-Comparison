function [x_GL_new, Ds_new] = rescale_diffMat_2D(x_GL, Ds, lb, ub)
% This function rescales the Gauss-Lobatto points and the differentiation
% matrices from [-1, 1] to [lb, ub].
%
% Args:
%       x_GL: Gauss-Lobatto points, each column corresponds to a component
%       Ds: Differentiation matices, a struct
%       lb: lower bounds
%       ub: upper bounds
% Return:
%       x_GL_new: Rescaled Gauss-Lobatto points
%       Ds_new: Rescaled differentiation matices

Ds_new = struct();

% x1
x1_GL = x_GL(:, 1);
D1 = Ds.D1;

x1_GL_new = (ub(1) - lb(1))/2 * (x1_GL + 1) + lb(1);
Ds_new.D1 = 2/(ub(1) - lb(1)) * D1;

% x2
x2_GL = x_GL(:, 2);
D2 = Ds.D2;

x2_GL_new = (ub(2) - lb(2))/2 * (x2_GL + 1) + lb(2);
Ds_new.D2 = 2/(ub(2) - lb(2)) * D2;

x_GL_new = [x1_GL_new, x2_GL_new];

