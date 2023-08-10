function A = rand_csymp(n,m)
%RANDCORTH  Random complex symplectic matrix.
%   2*N = matrix size, m = # of G-reflectors applied

I = eye(2*n); A = I;

for k = 1:m
    
    u = randn(2*n,1) + randn(2*n,1)*i;
    beta = randn + randn*i;
    J = [zeros(n) eye(n) ; -eye(n) zeros(n)];
    G = I + beta*u*u.'*J;
    A = A*G;

end