function f = CosineSquareDynamics(x,p)

f = p.k* cos(x).^2;
