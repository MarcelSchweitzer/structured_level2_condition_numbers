function A = construct_matrix(n,cnd,matrix,func)

% Randomly generate test matrix of chosen type. If func is logm or sqrtm,
% we make sure that the matrix does not have negative eigenvalues.

notgood = 1;

while notgood
    switch matrix      
        case  'orth'
            A = RandOrthMat(n);
        case 'symp'
            A = rand_rsymp(n/2,cnd);
        case 'perp'
            A = rand_rperp(n,cnd);
        case 'skew-sym'
            X = rand(n);
            A = (X-X')/2;
        case 'Hamilton'
            F1 = rand(n/2);
            G1 = rand(n/2);
            F = (F1+F1')/2;
            G = (G1+G1')/2;
            H = rand(n/2);
            
            A = [H, G; F, -H'];
        case 'triangular'
            D = diag(linspace(1e-2,1e-2*cnd,n));
            X = rand(n);
            A = X*D/X;
            A = -A;
            [~,T] = schur(A,'real');
            A = T;
    end
   
    if strcmp(func,'logm') == 1 || strcmp(func,'sqrtm') == 1      
        [Q_a,T_a] = schur(A,'real');
        
        % Check for any negative eigenvalues
        [~, T1] = rsf2csf(Q_a, T_a);
        
        notgood = any(imag(diag(T1)) == 0 & real(diag(T1)) <= 0 );
    elseif strcmp(func,'expm') == 1
        notgood = 0;
    end
    
end
