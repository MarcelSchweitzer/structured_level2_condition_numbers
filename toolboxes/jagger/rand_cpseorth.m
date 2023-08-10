function A = rand_cpseorth(p,q,m)
%RAND_CPSEORTH  Random complex pseudo-orthogonal matrix
%   m = # of G-reflectors applied

n = p + q; I = eye(n); A = I;

for k = 1:m
    
    u = randn(n,1) + randn(n,1)*i;
    sig = blkdiag(eye(p),-eye(q));
    H = I - (2/(u.'*sig*u))*u*u.'*sig;
    A = A*H;

end