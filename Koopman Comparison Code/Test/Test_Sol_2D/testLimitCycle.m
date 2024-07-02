%Test Result of Koopman Linearization for cosine square model 
% 
% Args:
%       dynamics: the dynamics f of the system, a function handle
%       p: inherent parameters of the system
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

clear;
clc;
%% Parameters
m = 200;
T = 20;
tspan = linspace(0, T, m+1);
n = 200;
N = 9;
x0 = [sqrt(2)/2, -sqrt(2)/2];
r = [sqrt(2)/8, sqrt(2)/8];
frac = 0.2;
op = 1;

p.k = 1;

fmodel = @(t, x, p) LimitCycle(t, x, p);
fdynamics = @(X, Y, p) LimitCycleDynamics(X, Y, p);
error=zeros(n+1);

%K
K = approximate_Koopman_2D(N,x_GL,Ds,dynamics,p);
eig(K)


