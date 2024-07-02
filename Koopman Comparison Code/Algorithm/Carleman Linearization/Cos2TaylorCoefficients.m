function cos2taylorCoefficients = Cos2TaylorCoefficients(x, N)
    cos2taylorCoefficients = zeros(1, N);
    
    for n = 0:N-1
        cos2taylorCoefficients(n+1) = (-1)^n * cos(x)^2 / factorial(2*n);
    end
end
