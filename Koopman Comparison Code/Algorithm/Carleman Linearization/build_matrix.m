function out = build_matrix(N,varargin)
    % build_matrix_N Construct a composite matrix with special structure
    %   Constructs a composite matrix with specific block-diagonal and
    %   off-diagonal structure involving Kronecker sums of input matrices F1 and F2.
    %   The structure of the composite matrix is governed by the input N.
    %
    %   Inputs:
    %       varargin - cofficient matrix F1,F2,...,Fk, and Fk is n by n^k
    %  
    %       N        - Integer value for the truncation order N
    %   Outputs:
    %       out      - Composite matrix with special structure

    %    The  carleman matrix will be constructed for each row, three parts
    %    are computed separatly left part, block part, right part.
    
    if nargin < 3
        error('at least 3 arguments required');
    end

    d = size(varargin{1}, 1); % d is the dimension of x

    % Initialize a cell array to store the output matrix blocks
    out = cell(1, N);

    % Iterate from 1 to N to build each block of the matrix
    for j = 1:N
        if nargin+j-2 > N             % last means the block part
            last = N-j+1;
        else                          
            last = nargin-1;
        end
        block_j = [];
        for i=1:last
            block_j = [block_j,kron_sum(varargin{i},j)];
        end
        left_j  = sparse(d^j, sum(arrayfun(@(i) d^i, 1:(j-1))));
        right_j = sparse(d^j, sum(arrayfun(@(i) d^i, (nargin+j-1):N)));
        
        out{j}   = [left_j, block_j, right_j];

        end
    
    % Concatenate the cell array vertically to form the final matrix
    out = vertcat(out{:});
    
end