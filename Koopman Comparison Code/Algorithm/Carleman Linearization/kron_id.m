function A = kron_id(n, k)
    % kron_id Compute the k-fold Kronecker product of the identity
    %   Computes the k-fold Kronecker product of the identity matrix of
    %   order n: I ⊗ I ⊗ ... ⊗ I, k times.
    %
    %   Inputs:
    %       n - integer representing the order (dimension of the identity)
    %       k - integer representing the power
    %   Outputs:
    %       A - Diagonal matrix with element 1 in the diagonal and order n^k

    % Create a diagonal matrix with ones and size n^k
    A = eye(n^k);
end


