function L = highFD(func, A, varargin)
%HIGHFD Computes high order Frechet derivatives of a matrix function.
% Calculates arbitrarily high order Frechet derivatives of a matrix
% function FUNC at the point A using [1, Alg. 4.2].
%
% Use the code as follows:
%   L = highFD(func, A, directions);
% where FUNC is a function handle to a matrix function,
% A is the matrix at which to evaluate the derivative and
% the directions Ei are placed in the last argument.
%
% For example using the matrix logarithm we could say
%   L = highFD(@logm, A, E1, E2, E3);
% to compute the third Frechet derivative in the directions E1, E2 and E3.
%
% References:
% [1] Nicholas J. Higham and Samuel D. Relton,
%     Higher Order Frechet Derivatives of Matrix Functions,
%     MIMS Eprint 2013.

E = varargin;
n = length(A);
k = length(E);

X = A;
for i = 1:k
    % kron is REALLY slow so trying this instead
    X = kronecker(eye(2), X) + kronecker( kronecker([0 1; 0 0], speye(2^(i-1))), E{i});
end
F = func(full(X));
L = F(1:n, end-n+1:end);

end