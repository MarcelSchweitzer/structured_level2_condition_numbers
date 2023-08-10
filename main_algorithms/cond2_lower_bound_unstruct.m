function [out] = cond2_lower_bound_unstruct(A, func, epsilon)
%COND2_LOWER_BOUND Lower bound on the level 2 condition number
%   A - Matrix at which to evaluate lvl2cond
%   func - Function which calculates f(A)

n = length(A);
maxiter = 3000;

% Compute level-1 condition number
condfA = cond_level1(A, func);

% Loop through matrix types to ensure Z is a structured perturbation
cond2fA_obj = @(z) -abs(cond_level1(A + epsilon*((reshape(z,n,n))/norm(reshape(z,n,n), 'fro')),func) - condfA)/epsilon;

z0 = randn(n^2, 1);
[~, out] = fminsearch(cond2fA_obj, z0, optimset('MaxFunEvals', maxiter, 'MaxIter', maxiter, 'TolX', 0, 'TolFun',0));%,'PlotFcns',@optimplotfval));
out = -1 * out;

% SUBFUNCTIONS
    function [O] = cond_level1(X, func)
        K = highKronForm(func, X, 1);

        % Compute level-1 condition number
        O = norm(K, 2);
    end


end