function B = qmult_unit(A,method)
%QMULT Pre-multiply matrix by random unitary matrix.
%   QMULT_UNIT(A) returns Q*A where Q is a random unitary matrix
%   from the Haar distribution of dimension the number of rows in A.
%   Special case: if A is a scalar then QMULT(A) is the same as QMULT(EYE(A)).
%   QMULT_UNIT(A,METHOD) specifies how the computations are carried out.
%   METHOD = 0 is the default, while METHOD = 1 uses a call to QR, which
%   is much faster for large dimensions, even though it uses more flops.

[n, m] = size(A);

%  Handle scalar A.
if max(n,m) == 1
   n = A;
   A = eye(n);
end

if nargin == 2 & method == 1
   F = randn(n) + randn(n)*i;
   [Q,R] = qr(F);
   B = Q*diag(exp(rand(n,1)*2*pi*i))*A;
   return
end

for k = n-1:-1:1
    
    % Generate random complex unitary Hermitian Householder.
    x = randn(n-k+1,1) + randn(n-k+1,1)*i;
    [v, beta] = gallery('house',x);

    % Apply the transformation to A.
    y = v'*A(k:n,:);
    A(k:n,:) = A(k:n,:) - beta*v*y;
    
end

A = diag(exp(rand(n,1)*2*pi*i))*A;
B = A;