function [x, D] = chebyshevDiff(N, a, b)
% This function computes the Gauss-Chebyshev-Lobatto 
% differentiation matrix in the interval [a,b]
%
% Input N   -> polynomial order 
%       a,b -> interval [a,b]
%
% Output x  -> vector of (N+1) Chebyshev points  in [a,b].
%        D  -> (N+1) x (N+1) diff matrix such that df/dx = D*f in [a,b]
%

if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)'; 
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';                  
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries

x = flipud(x);
D = -D;

D = 2/(b-a)*D;
x = (b-a)/2*(x+1) + a;
