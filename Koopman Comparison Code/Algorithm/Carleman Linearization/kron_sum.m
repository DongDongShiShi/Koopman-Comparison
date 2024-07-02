function A = kron_sum(F, k)
    % kron_sum Compute the Kronecker sum of order k
    %   Computes the Kronecker sum of order k as described in the input documentation.
    %   The function computes the sum as F ⊕ F ⊕ ... ⊕ F with a total of k summands.
    %
    %   Inputs:
    %       F - matrix
    %       k - integer representing the order of the sum
    %   Outputs:
    %       A - The Kronecker sum

    if k < 1
        error('expected k ≥ 1, got %d', k);
    elseif k == 1
        A = F;
        return;
    end

    n = size(F, 1); %  the first dimension of the matrix 
    k1 = 0; % the number of the terms on the left of F
    k2 = k - 1; % the number of the terms on the right of F
    A = kron_sandwich(F, n, k1, k2);  % I^(⊗k1) ⊗ F ⊗ I^(⊗k2) the current term computed
    
    % compute the sum
    for i = 2:k
        k1 = k1 + 1;
        k2 = k2 - 1;
        B = kron_sandwich(F, n, k1, k2);
        A = A + B;
    end
end
