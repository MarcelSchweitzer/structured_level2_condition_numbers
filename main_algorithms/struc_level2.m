function [x,r,l_2,sl_2,slb_2,lb_2, CN]= struc_level2(n,matrix,f,func,strct,group)

% function STRUC_LEVEL2 compares the structured/unstructured level-2
% condition number for some specific matrix functions.

% Matrix:
%      Orthogonal       'orth'
%      Symplectic       'symp'
%      Perplectic       'perp'
%      Skew-symmetric   'skew-sym'
%      Hamiltonian      'Hamilton'
%
% f: Matrix square root (@sqrtm), matrix exponential (@expm), matrix logarithm (@logm)
% func: 'sqrtm', 'expm', 'logm'
% STRCT gives the choice of the matrix M:
%      'I', 'J', 'R', 'Sigma'
% GROUP specifies:
%      Automorphism group ('auto'), Jordan ('jord') and Lie algebra ('lie')


CN = [2 5 1e1 5e1 1e2 5e2 1e3 5e3 1e4 5e4 1e5 5e5 1e6 5e6 1e7 1e8  1e9  1e10];
num_samples = 2; % number of test matrices to average over for each data point

condA = zeros(length(CN),1);
l_2_arr = zeros(num_samples,length(CN));
sl_2_arr = zeros(num_samples,length(CN));
r_arr = zeros(num_samples,length(CN));
slb_2_arr = zeros(num_samples,length(CN));
lb_2_arr  = zeros(num_samples,length(CN));

for k = 1:length(CN)
    disp(k)
    for s = 1:num_samples
        cnd = CN(k);

        % Build test matrices.
        A = construct_matrix(n,cnd,matrix,func);

        N = n^2;
        I_N = eye(N);

        condA(k) = cond(A);

        % Compute the second Kronecker form of FUNC at A.
        K = highKronForm(f, A, 2);

        % Construct the basis of the tangent spaces.
        Q = construct_base(A,strct,group);

        % Upper bounds of the unstructured/structured level-2 condition number.
        l_2_arr(s,k) = norm(K);
        sl_2_arr(s,k) =  norm(kron((Q*pinv(Q))',I_N)*K*(Q*pinv(Q)));

        r_arr(s,k) = sl_2_arr(s,k)/l_2_arr(s,k);

        % Lower bound on structured level-2 condition number.
        slb_2_arr(s,k) = cond2_lower_bound(A, f, strct, group, matrix);

        % Lower bound on unstructured level-2 condition number.
        if strcmp(matrix,'triangular') == 1
            epsilon = 1e-9;
        else
            epsilon = 1e-12;
        end
        lb_2_arr(s,k) = cond2_lower_bound_unstruct(A, f, epsilon);
    end
end

if num_samples > 1
    l_2 = sum(l_2_arr)/num_samples;
    sl_2 = sum(sl_2_arr)/num_samples;
    r = sum(r_arr)/num_samples;
    lb_2 = sum(lb_2_arr)/num_samples;
    slb_2 = sum(slb_2_arr)/num_samples;
else
    l_2 = l_2_arr;
    sl_2 = sl_2_arr;
    r = r_arr;
    lb_2 = lb_2_arr;
    slb_2 = slb_2_arr;
end
x = 1:length(CN);