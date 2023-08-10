function A = rand_pseunit(p,q,c,symm,method)
%RAND_PSEUNIT  Random complex pseudo-unitary matrix.
%   A = RAND_PSEUNIT(P,Q,C) forms a random (P+Q)-by-(P+Q) 
%   complex pseudo-unitary matrix A, with COND(A) = C.
%   Pseudo-unitary means that A'*SIGMA*A = SIGMA, with A complex,
%   where SIGMA = BLKDIAG(EYE(P),-EYE(Q)).
%   If omitted, C defaults to SQRT(1/EPS).
%
%   RAND_PSEUNIT(N) and RAND_PSEUNIT(N,[],C)
%   both produce an N-by-N matrix with P = CEIL(N/2), Q = FLOOR(N/2).
%
%   The full calling sequence is RAND_PSEUNIT(P,Q,C,SYMM,METHOD).
%   If SYMM is nonzero symmetry is enforced.
%   The argument METHOD specifies how the underlying unitary
%   transformations are carried out.  If METHOD is nonzero
%   a call to QR is used, which is much faster than the default
%   method for large dimensions, though it uses more flops.
%   August 15, 17:30

if nargin < 2 || isempty(q), q = floor(p/2); p = p-q; end
if nargin < 3 || isempty(c), c = sqrt(1/eps); end
if nargin < 4 || isempty(symm), symm = 0; end
if nargin < 5, method = 0; end

% This function requires q >= p, so...
if p > q
   A = randjunit(q,p,c,symm,method); % diag(eye(q),-eye(p))-unitary matrix.
   A = A([q+1:p+q 1:q],:); % Permute to produce pseudo-unitary matrix.
   A = A(:,[q+1:p+q 1:q]);
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
   U = blkdiag(qmult_unit(p, method),qmult_unit(q, method));
   A = U*A*U';
   A = (A + A')/2; % Ensure matrix is symmetric.
   return
end

A = left_mult(A,p,q,method); % Left multiplications by unitary matrices.
A = left_mult(A',p,q,method); % Right multiplications by unitary matrices.

function A = left_mult(A,p,q,method)
%LEFT_MULT   Left multiplications by random unitary matrices.
A(1:p,:) = qmult_unit(A(1:p,:), method);
A(p+1:p+q,:) = qmult_unit(A(p+1:p+q,:), method);