function U = construct_base(X,strct,group)

% CONSTRUCT_BASE function gives the base of the tangent space
% of 'group' at X.
% STRCT gives the choice of the matrix.
% GROUP specifies:
% automorphism group (G), Jordan (J) and Lie algebra (L)

n = length(X);
nn = n^2;
I = eye(n);

switch strct
    case 'I'
        mu = 1;
        M = eye(n);
    case 'J'
        mu = -1;
        M = [zeros(n/2) eye(n/2); -eye(n/2) zeros(n/2)];
    case 'R'
        mu = 1;
        M = zeros(n); M(n:n-1:n^2-n+1) = 1;
    case 'S1'
        mu = 1;
        p = n/2; q = n-p;
        M = blkdiag(eye(p),-eye(q));
    case 'S2'
        mu = 1;
        p = 1; q = n-p;
        M = blkdiag(eye(p),-eye(q));
    case 'T'
        % triangular matrix
        d = (diag(X,-1) ~= 0);
        nd = n*(n+1)/2 + nnz(d);
        U = zeros(n^2,nd);
        s = 0;
        e = zeros(n,1);
        
        for j = 1:n
            if j < n && d(j) == 1
                s = s+1;
                ej  = e; ej(j)    = 1;
                ej1 = e; ej1(j+1) = 1;
                U(:,s) = reshape(ej1*ej',n^2,1);
            end
            ej = e; ej(j) = 1;
            for i = 1:j
                s = s+1;
                ei  = e; ei(i) = 1;
                U(:,s) = reshape(ei*ej', n^2, 1);
            end
        end
        return %nothing more to do here, skip rest of function
end

sL = n*(n-mu)/2;
sJ = n*(n+mu)/2;

DL = zeros(nn,sL);
DJ = zeros(nn,sJ);

k = 1;
for j = 1:n
    for i = 1:j-1
        DL((i-1)*n+j,k) = 1/sqrt(2);
        DL((j-1)*n+i,k) = -mu/sqrt(2);
        DJ((i-1)*n+j,k) = 1/sqrt(2);
        DJ((j-1)*n+i,k) = mu/sqrt(2);
        k = k+1;
    end
end

if mu == -1
    for i=1:n
        DL((i-1)*n+i,k) = 1;
        k = k+1;
    end
else
    for i=1:n
        DJ((i-1)*n+i,k) = 1;
        k = k+1;
    end
end

switch group
    case 'auto'
        XX = mu*X*M;
        U = kron(I,XX)*DL;
    case 'jord'
        U = mu*kron(I,M)*DJ;
    case 'lie'
        U = mu*kron(I,M)*DL;
end

