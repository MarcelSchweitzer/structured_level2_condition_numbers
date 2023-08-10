function K = highKronForm(func, A, k)
% KRONFORM Computes the kth Kronecker form of FUNC at A.
% Calculates arbitrarily high Kronecker forms for arbitrarily high
% Frechet derivatives using [1, Alg. 5.2].
%
% This code is called as:
%   K = highKronForm(func, A, k);
% where FUNC is a function handle to a matrix function such as expm,
% A is the matrix at which to evaluate the Kronecker form and
% k is the order of the Kronecker form to be returned.
%
% To see the relationship between the Frechet derivative and the
% Kronecker form see [1, (5.3)].
%
% References:
% [1] Nicholas J. Higham and Samuel D. Relton,
%     Higher Order Frechet Derivatives of Matrix Functions,
%     MIMS Eprint 2013.

n = length(A);
recursion_level = 0;
E = cell(1,k); % Stores the k direction matrices

K = kronrecursion;

% Nested function

    function M = kronrecursion()
    % Recursive until the kth level when the FD is calculated
    
        % Check if we are on the last level
        if recursion_level == k
            % Return the derivative
            M = highFD(func, A, E{:});
        else
            % Continue with recursion
            M = zeros(n^(2*(k-recursion_level)),n^2);
            recursion_level = recursion_level + 1; % Current depth
            for j = 1:n^2
                % Create direction matrix
                e = zeros(n^2,1); e(j) = 1;
                myE = zeros(n); myE(:) = e;
                E{recursion_level} = myE;
                % Go to next level of recursion
                L = kronrecursion;
                M(:,j) = L(:);
            end
            recursion_level = recursion_level - 1; %Finished this level
        end
    end

end