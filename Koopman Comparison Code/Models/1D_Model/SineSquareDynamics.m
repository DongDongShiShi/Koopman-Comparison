function f = SineSquareDynamics(x,p)

f = p.k* sin(x).^2;
