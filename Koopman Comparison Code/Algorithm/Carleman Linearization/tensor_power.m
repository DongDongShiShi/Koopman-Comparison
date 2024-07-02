function result = tensor_power(x0, N)
% Input vector x0
% Input integer N

%Output vector x0 tensor power N

result = x0;
    for k = 2:N
        result = kron(result, x0);
    end
end