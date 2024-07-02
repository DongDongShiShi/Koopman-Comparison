function [init1,init2] = compute_initial_2D(x_GL)
% 1inear combination of eigenmatrices. 
% This linear combination projects the obvervable into the vector space 
% spanned by eigenmatrices.
% 
% Args:
%       x_GL: Gauss-Lobatto points
% Returns:
%       init1,init2: Initial Condition with observable as g(x)=x

[m, ~] = size(x_GL);

x1_GL = x_GL(:, 1);
x2_GL = x_GL(:, 2);

% c1 = V \ repmat(x1_GL, m, 1);
% c2 = V \ repelem(x2_GL, m);
init1 = repmat(x1_GL, m, 1);
init2 = repelem(x2_GL, m);


