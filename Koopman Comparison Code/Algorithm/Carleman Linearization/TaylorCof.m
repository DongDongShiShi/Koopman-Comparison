% Polynomial Cofficient by Taylor Expansion
function out = TaylorCof(f, x0, k)
% This function is used to compute cofficient of Taylor series, notice we
% obtain simplified result. We transfer taylor series at x0 as
% F_0 + F_1*x + F_2*x^2 + F_3*x^3 + ... + F_k*x^k
% 
% math formula based on Taylor Series is that 
% F_j= sum(i = j to k) (f^(i)(x0)/i!)* C_i^{i-j}*(-x0)^{i-j}
% -------------------------------------------------------------------------------
% Input ---- f : 1D function handle
%       ---- x0: real value, the working point 
%       ---- k : integer, the expansion order
%
% Output----out: A vector to store all the cofficient of monomials
%                out=[F_0,F_1,F_2,...,F_k]
%




    out = cell(1,k);  % Initialize a cell array to store k matrix expressions
    syms x;  % Declare the symbolic variable x
    diff_term=zeros(1,k);
    for n=1:k
     diff_term(n)= double(subs(diff(f,x,n),x,x0));   %Compute n-th order differentiation f^(n)(x0) and store it for future use
    end
    for j = 1:k
        out{j} = 0;  % Initialize the j-th matrix expression to 0
        
        for i = j:k
            
            % x_term (-x0)^{i-j}  to aviod 0^0=NaN
            if x0==0 && i-j==0
                x_term=1;
            else
                x_term = (-x0)^(i - j);
            end
            % Update the j-th matrix expression 
            %out{j} = out{j} + double(diff_term(i)*(factorial(i))^2 / (factorial(i - j) * factorial(j)) * x_term);
            out{j} = out{j} + diff_term(i)*(factorial(i))^2 / (factorial(i - j) * factorial(j)) * x_term;
        end
    end
end
