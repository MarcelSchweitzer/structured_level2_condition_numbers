function A = rand_rsymp(n,c,method)
%RAND_RSYMP  Random real symplectic matrix.
%   A = RAND_RSYMP(N,C,METHOD) forms a random 2N-by-2N real symplectic 
%   matrix A, where C determines COND(A) via three possible cases:
%    - If C is a scalar, then COND(A) = C.
%    - If C is a real N-vector with all elements >= 1, C specifies 
%      SIGMA(1) >= SIGMA(2) >= ... >= SIGMA(N), the first N singular
%      values of A. The remainding singular values are calculated 
%      according to constraints: SIGMA(N+K) = 1/SIGMA(K).
%    - If omitted, C defaults to SQRT(1/EPS).
%   Real symplectic means that A'*J*A = J, where 
%   J = [ZEROS(N),EYE(N);-EYE(N),ZEROS(N)].
%
%   The argument METHOD specifies how the random complex unitary
%   matrices are generated. If METHOD is nonzero, a call to QR is used, 
%   which is much faster than the default method for large dimensions, 
%   though it uses more flops.
%   Aug 17, 16:30

if nargin < 2 || isempty(c), c = sqrt(1/eps); end
if nargin < 3, method = 0; end

s = zeros(1,2*n);
p = length(c);

switch p
    
    case 1
    if c >= 1
        s(1) = sqrt(c);
        s(n:-1:2) = sort(1 + (s(1)-1)*rand(n-1,1)); 
        s(n+1:2*n) = 1./s(1:n);
    else
        error('Illegal value for C.  To specify COND set C >= 1.')
    end
    
    case n
    c = sort(c);
    if c(1) < 1
        error('Illegal value for C. All values must be >=1.')
    else
        s(1:n) = c(n:-1:1);
        s(n+1:2*n) = 1./s(1:n);
    end
    
    otherwise
    error('Vector of singular values must be of length 1 or N.')
    
end   

U = symp_orth(n,method); V = symp_orth(n,method);
A = U*diag(s)*V';

% subfunction
function U = symp_orth(n,method)
%SYMP_ORTH  Random symplectic orthogonal matrix
B = qmult_unit(n,method);
Re = real(B); Im = imag(B);
U = [Re,Im;-Im,Re];