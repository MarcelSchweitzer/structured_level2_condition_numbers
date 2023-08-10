function A = rand_rperp(n,c,method)
%RAND_RPERP  Random real perplectic matrix.
%   A = RAND_RPERP(N,C,METHOD) forms a random N-by-N real perplectic 
%   matrix A, where COND(A) = C.
%   If omitted, C defaults to SQRT(1/EPS).
%   Real perplectic means that A'*R*A = R, where R defined by
%   R = ZEROS(N); R(N:N-1:N^2-N+1) = 1 (R = 1 if N = 1)
%   The argument METHOD specifies how the random orthogonal
%   matrices are generated. If METHOD is nonzero, a call to QR is used, 
%   which is much faster than the default method for large dimensions, 
%   though it uses more flops.
%   August 29, 13:30

if nargin < 2 || isempty(c), c = sqrt(1/eps); end
if nargin < 3, method = 0; end

p = floor(n/2); s = zeros(1,n);

if c < 1 
   error('Illegal value for C.  To specify COND set C >= 1.')
   
elseif rem(n,2)
   
   s(1) = sqrt(c);
   s(p:-1:2) = sort(1 + (s(1)-1)*rand(p-1,1));
   s(p+1) = mysign(randn);
   s(p+2:n) = 1./s(p:-1:1);

   U = perp_orth_odd(p,method); V = perp_orth_odd(p,method);
   A = U*diag(s)*V';

else
    
   s(1) = sqrt(c);
   s(p:-1:2) = sort(1 + (s(1)-1)*rand(p-1,1));
   s(p+1:n) = 1./s(p:-1:1);

   U = perp_orth_even(p,method); V = perp_orth_even(p,method);
   A = U*diag(s)*V';
   
end

% subfunction 1
function U = perp_orth_even(k,method)
%PERP_ORTH_2N   Random perplectic orthogonal of size 2k
P = gallery('qmult',k,method); Q = gallery('qmult',k,method);
X = P + Q(k:-1:1,k:-1:1); Y = P(:,k:-1:1) - Q(k:-1:1,:);
U = (1/2)*[ X  Y ; Y(k:-1:1,k:-1:1) X(k:-1:1,k:-1:1) ];

% subfunction 2
function U = perp_orth_odd(k,method)
%PERP_ORTH_ODD2  Random perplectic orthogonal of size 2k+1
P = gallery('qmult',k+1,method); Q = gallery('qmult',k,method);
U = zeros(2*k+1);
U = (1/2)*[ P(1:k,1:k) + Q(k:-1:1,k:-1:1) sqrt(2)*P(1:k,k+1) ...
            P(1:k,k:-1:1) - Q(k:-1:1,:) ; 
            sqrt(2)*P(k+1,1:k) 2*P(k+1,k+1) sqrt(2)*P(k+1,k:-1:1) ; ...
            P(k:-1:1,1:k) - Q(:,k:-1:1) sqrt(2)*P(k:-1:1,k+1) ...
            P(k:-1:1,k:-1:1) + Q ];