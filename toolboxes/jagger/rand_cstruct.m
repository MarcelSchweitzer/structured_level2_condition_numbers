function A = rand_cstruct(f,n,c)
%RAND_CSTRUCT   Random complex structured matrices
%  * A = RAND_CSTRUCT(1,N,C) forms a random N-by-N complex orthogonal 
%    matrix A, with COND(A) approximately equal to C.
%  * A = RAND_CSTRUCT(2,[P Q],C) forms a random (P+Q)-by-(P+Q) 
%    complex pseudo-orthogonal matrix A, with COND(A) approximately 
%    equal to C.
%  * A = RAND_CSTRUCT(3,N,C) forms a random N-by-N complex symplectic
%    matrix A, with COND(A) approximately equal to C.
%  * A = RAND_CSTRUCT(4,N,C) forms a random N-by-N conjugate symplectic
%    matrix A, with COND(A) approximately equal to C.
%       
%    In each case, X contains the LS solution [A0,A1,A2,A3] for
%    LOG(C) = A0 + A1*K + A2*N + A3*K^2 from empirical data, where
%    K is the number of G-reflectors applied, N is the matrix size, and
%    C is the observed condition number of the matrix.
%
%    QUAD contains the co-efficients of the quadratic 
%    A3*K^2 + A1*K + (A2*N + A0 - LOG(C)). 
%    QUAD is solved for K, and K is rounded to the nearest integer to
%    formulate how many G-reflectors to apply.

if nargin < 3 || isempty(c), c = sqrt(1/eps); end

switch f
    
    case 1
        
    I = eye(n); A = I;
    X = [ 2.0344 0.0044 1.9239 -0.0249 ];
    quad = [ X(4), X(3), X(2)*n + X(1) - log(c) ];
    k = round(min(sort(roots(quad))));
    if k < 1, k = 1; end
    fprintf('Number of G-reflectors being applied : %g', k)
    
    for j = 1:k
    
        u = randn(n,1) + randn(n,1)*i;
        G = I - (2/(u.'*u))*u*u.';
        A = A*G;
    
    end
    
    case 2

    [a,b] = size(n);
    if a*b ~= 2
        error('Must enter vector [p q] for matrix size.')
    end
    
    p = n(1); q = n(2); n = p + q; I = eye(n); A = I;
    X = [ 1.9510 0.0058 1.9080 -0.0245 ]; 
    quad = [ X(4), X(3), X(2)*n + X(1) - log(c) ]; 
    k = round(min(sort(roots(quad))));
    if k < 1, k = 1; end
    fprintf('Number of G-reflectors being applied : %g', k)
    
    for j = 1:k
    
        u = randn(n,1) + randn(n,1)*i;
        sig = blkdiag(eye(p),-eye(q));
        G = I - (2/(u.'*sig*u))*u*u.'*sig;
        A = A*G;

    end
    
    case 3

    I = eye(2*n); A = I;
    X = [ 3.9794 0.0249 5.8397 -0.2311 ];
    quad = [ X(4), X(3), X(2)*n + X(1) - log(c) ];
    k = round(min(sort(roots(quad))));
    if k < 1, k = 1; end
    fprintf('Number of G-reflectors being applied : %g', k)
    
    for j = 1:k
    
        u = randn(2*n,1) + randn(2*n,1)*i;
        beta = randn + randn*i;
        J = [zeros(n) eye(n) ; -eye(n) zeros(n)];
        G = I + beta*u*u.'*J;
        A = A*G;
    
    end
    
    case 4
    
    I = eye(2*n); A = I;
    X = [ 3.6339 0.0063 2.0899 -0.0274]; 
    quad = [ X(4), X(3), X(2)*n + X(1) - log(c) ];
    k = round(min(sort(roots(quad))));
    if k < 1, k = 1; end
    fprintf('Number of G-reflectors being applied : %g', k)
    
    for j = 1:k
    
        u = randn(2*n,1) + randn(2*n,1)*i;
        J = [zeros(n) eye(n) ; -eye(n) zeros(n)];
        q = u'*J*u; r = -1/q; r = r - real(r);
        beta = r*exp(2*pi*i*rand) + r;
        G = I + beta*u*u'*J;
        A = A*G;

    end
    
    otherwise
    error('Must choose matrix type 1-4.')
        
end