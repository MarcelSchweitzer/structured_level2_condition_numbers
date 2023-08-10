function A = rand_rperp2(n,c,method)
%RAND_RPERP  Random real perplectic matrix.
%   A = RAND_RPERP(N,C,METHOD) forms a random N-by-N real perplectic 
%   matrix A, where COND(A) = C.
%   If omitted, C defaults to SQRT(1/EPS).
%   Real perplectic means that A'*R*A = R, where R defined by
%   R = ZEROS(N,N); R(N:N-1:N^2-N+1) = 1 (R = 1 if N = 1)
%   The argument METHOD specifies how the random orthogonal
%   matrices are generated. If METHOD is nonzero, a call to QR is used, 
%   which is much faster than the default method for large dimensions, 
%   though it uses more flops.
%   August 29, 13:30
%   need to edit comments if use it - maybe sort into the entire even case
%   and the entire odd case - input vector of p+1 for odd, with (p+1)-th
%   being a plus or minus 1?

if nargin < 2 || isempty(c), c = sqrt(1/eps); end
if nargin < 3, method = 0; end

p = floor(n/2); s = zeros(1,n);
f = length(c);

switch f
    
    case 1
    if c < 1 
        error('Illegal value for C.  To specify COND set C >= 1.')
    elseif rem(n,2)
        s(1) = sqrt(c);
        s(p:-1:2) = sort(1 + (s(1)-1)*rand(p-1,1));
        s(p+1) = mysign(randn);
        s(p+2:n) = 1./s(p:-1:1);
    else
        s(1) = sqrt(c);
        s(p:-1:2) = sort(1 + (s(1)-1)*rand(p-1,1));
        s(p+1:n) = 1./s(p:-1:1);
    end

    case p
    c = sort(c);
    if c(1) < 1
        error('Illegal value for C. All values must be >=1.')
    elseif rem(n,2)
        s(1:p) = c(p:-1:1);
        s(p+1) = mysign(randn); % this is still random hmmm
        s(p+2:n) = 1./s(p:-1:1);
    else
        s(1:p) = c(p:-1:1);
        s(p+1:n) = 1./s(p:-1:1);
    end
    
    otherwise
    error('Vector of singular values must be of length 1 or FLOOR(N/2).')
    
end     

if rem(n,2)
    
    U = perp_orth_odd(p,method); V = perp_orth_odd(p,method);
    A = U*diag(s)*V';

else

    U = perp_orth_even(p,method); V = perp_orth_even(p,method);
    A = U*diag(s)*V';
    
end

% subfunction 1
function U = perp_orth_even(k,method)
%PERP_ORTH_2N   Random perplectic orthogonal of size 2k
P = gallery('qmult',k,method); Q = gallery('qmult',k,method);
if k == 1, R = 1; else R = zeros(k); R(k:k-1:k^2-k+1) = 1; end
X = P + R*Q*R; Y = P*R - R*Q;
U = (1/2)*[ X  Y ; R*Y*R R*X*R ];

% subfunction 2
function U = perp_orth_odd(k,method)
%PERP_ORTH_ODD2  Random perplectic orthogonal of size 2k+1
P = gallery('qmult',k+1,method); Q = gallery('qmult',k,method);
U = zeros(2*k+1);
U = (1/2)*[ P(1:k,1:k) + Q(k:-1:1,k:-1:1) sqrt(2)*P(1:k,k+1) ...
            P(1:k,k:-1:1) - Q(k:-1:1,1:k) ; 
            sqrt(2)*P(k+1,1:k) 2*P(k+1,k+1) sqrt(2)*P(k+1,k:-1:1) ; ...
            P(k:-1:1,1:k) - Q(1:k,k:-1:1) sqrt(2)*P(k:-1:1,k+1) ...
            P(k:-1:1,k:-1:1) + Q ];