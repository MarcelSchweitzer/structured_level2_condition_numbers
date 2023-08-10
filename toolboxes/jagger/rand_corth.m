function A = rand_corth(n,m)
%RANDCORTH  Random complex orthogonal matrix.
%   N = matrix size, m = # of G-reflectors applied

I = eye(n); A = I;

for k = 1:m
    
    u = randn(n,1) + randn(n,1)*i;
    H = I - (2/(u.'*u))*u*u.';
    A = A*H;
    
end