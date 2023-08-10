function [out] = cond2_lower_bound(A, func, strct, group, matrix)
%COND2_LOWER_BOUND Lower bound on the level 2 condition number
%   A - Matrix at which to evaluate lvl2cond
%   func - Function which calculates f(A)
%   strct - Choice of M matrix
%   group - Which group structure

B = construct_base(A,strct,group);

n = length(A);
maxiter = 3000;

condfA = cond_level1(A, func, strct, group);

if strcmp(matrix,'triangular')
    epsilon = 1e-9;
else
    epsilon = 1e-3;
end

[~,b] = size(B);
z0 = randn(b, 1);
cond2fA_obj = @(z) -abs(cond_level1(A + epsilon*reshape((B*z)/norm(B*z),n,n), func, strct, group) - condfA)/epsilon;

[~, out] = fminsearch(cond2fA_obj, z0, optimset('MaxFunEvals', maxiter, 'MaxIter', maxiter, 'TolX', 0, 'TolFun',0));
out = -1 * out;

% SUBFUNCTIONS
    function [O] = cond_level1(X, func, strct, group)
        % Basis of structure tangent space
        Ml = construct_base(X, strct, group);
        Kl = highKronForm(func, X, 1);

        % Compute level-1 condition number
        O = norm(Kl*Ml*pinv(Ml), 2);
    end

end