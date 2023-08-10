function [x,l_2,sl_2,slb_2,lb_2, condA]= struc_level2_benchmark()

% function STRUC_LEVEL2_BENCHMARK compares the structured/unstructured level-2
% condition number for the matrix exponential, using several test matrices
% from benchmark collections

f = @expm;
load 'data/triangular_benchmark_matrices.mat' triangular_benchmark_matrices

l_2 = zeros(length(triangular_benchmark_matrices),1);
sl_2  = zeros(length(triangular_benchmark_matrices),1);
r = zeros(length(triangular_benchmark_matrices),1);
slb_2  = zeros(length(triangular_benchmark_matrices),1);
lb_2  = zeros(length(triangular_benchmark_matrices),1);
condA = zeros(length(triangular_benchmark_matrices),1);

for k = 1:length(triangular_benchmark_matrices)
    disp(k)

    % Load test matrices.
    A = triangular_benchmark_matrices{k};
    n = size(A,1);
    N = n^2;
    I_N = eye(N);

    condA(k) = cond(A);

    % Compute the second Kronecker form of FUNC at A.
    K = highKronForm(f, A, 2);

    % Construct the basis of the tangent spaces.
    Q = construct_base(A,'T', []);

    % Upper bounds of the unstructured/structured level-2 condition number.
    l_2(k) = norm(K);
    sl_2(k) =  norm(kronecker((Q*pinv(Q))',I_N)*K*(Q*pinv(Q)));

    % Lower bound on structured level-2 condition number.
    slb_2(k) = cond2_lower_bound(A, @expm, 'T', [], 'triangular');

    % Lower bound on unstructured level-2 condition number.
    lb_2(k) = cond2_lower_bound_unstruct(A, @expm, 1e-9);
end

x = 1:length(triangular_benchmark_matrices);