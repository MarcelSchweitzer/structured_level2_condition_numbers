function A = randjorth(p,q,c,symm,method)
%RANDJORTH  Random J-orthogonal matrix.
%   A = GALLERY('RANDJORTH',P,Q,C) forms a random (P+Q)-by-(P+Q) J-orthogonal
%   matrix A, where J = BLKDIAG(EYE(P),-EYE(Q)) and COND(A) = C.
%   J-orthogonality means that A'*J*A = J, and such matrices are
%   sometimes called hyperbolic.
%   If omitted, C defaults to SQRT(1/EPS).
%
%   GALLERY('RANDJORTH',N) and GALLERY('RANDJORTH',N,[],C)
%   both produce an N-by-N matrix with P = CEIL(N/2), Q = FLOOR(N/2).
%
%   The full calling sequence is
%        GALLERY('RANDJORTH',P,Q,C,SYMM,METHOD).
%   If SYMM is nonzero symmetry is enforced and a symmetric positive
%   definite matrix is produced.
%   The argument METHOD specifies how the underlying orthogonal
%   transformations are carried out.   If METHOD is nonzero
%   a call to QR is used, which is much faster than the default
%   method for large dimensions, though it uses more flops.

%   Reference:
%   N. J. Higham, J-orthogonal matrices: Properties and generation,
%   Numerical Analysis Report No. 408, Manchester Centre for
%   Computational Mathematics, Manchester, England, Sept. 2002.
%
%   Nicholas J. Higham

if nargin < 2 || isempty(q), q = floor(p/2); p = p-q; end
if nargin < 3 || isempty(c), c = sqrt(1/eps); end
if nargin < 4 || isempty(symm), symm = 0; end
if nargin < 5, method = 0; end

% This function requires q >= p, so...
if p > q
   A = randjorth(q,p,c,symm,method); % diag(eye(q),-eye(p))-orthogonal matrix.
   A = A( [q+1:p+q 1:q], :); % Permute to produce J-orthogonal matrix.
   A = A(:, [q+1:p+q 1:q]);
   return
end

if c >= 1
   c(1) = (1+c)/(2*sqrt(c));
   c(2:p) = 1 + (c(1)-1)*rand(p-1,1);
elseif p ~= 1
   error('Illegal value for C.  To specify COND set C >= 1.')
end

s = sqrt(c.^2-1);

A = blkdiag([diag(c) -diag(s); -diag(s) diag(c)], eye(q-p));

if symm
   U = blkdiag(gallery('qmult',p, method),gallery('qmult',q, method));
   A = U*A*U';
   A = (A + A')/2;       % Ensure matrix is symmetric.
   return
end

A = left_mult(A,p,q,method);  % Left multiplications by orthogonal matrices.
A = left_mult(A',p,q,method); % Right multiplications by orthogonal matrices.

function A = left_mult(A,p,q,method)
%LEFT_MULT   Left multiplications by random orthogonal matrices.
A(1:p,:) = gallery('qmult',A(1:p,:), method);
A(p+1:p+q,:) = gallery('qmult',A(p+1:p+q,:), method);