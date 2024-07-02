function A = kron_sandwich(F, n, k1, k2)
    % kron_sandwich Compute A ⊗ F ⊗ C where A = I^{⊗ k1} and C = I^{⊗ k2}
    %   Computes the Kronecker product where A and C are k1-fold and k2-fold 
    %   Kronecker products of the identity of order n, and F is in between.
    %
    %   Inputs:
    %       F  - matrix
    %       n  - integer, dimension of the identity
    %       k1 - nonnegative integer
    %       k2 - nonnegative integer
    %      ！！！ F could not be a square matrix, but for its size n*m.  
    %      ！！！ where n should be exactly the same as the input n

    %   Outputs:
    %       A - The Kronecker product I^{⊗ k1} ⊗ F ⊗ I^{⊗ k2}

    % Compute A ⊗ F
    if k1 == 0
        AF = F;
    else
        A = kron_id(n, k1);
        AF = kron(A, F);
    end

    % Compute A ⊗ F ⊗ C
    if k2 == 0
        A = AF;
    else
        C = kron_id(n, k2);
        A = kron(AF, C);
    end
end
