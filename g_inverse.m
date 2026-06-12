function [x, rmse] = g_inverse(A, b, avoid_neg, accy, quad_method)
% g_inverse solves the tomographic linear system and returns fit RMSE.
% It uses constrained quadratic optimization, NNLS, or lsqminnorm depending
% on the selected options. MATLAB lsqnonneg is avoided because it is too slow here.

nVoxels = size(A,2);
At = transpose(A);

if quad_method > 0
    
    % Constrained quadratic optimization using MATLAB quadprog.
    N = At * A;
    f = -At * b;
    N=(N+N')/2;

    % interior-point-convex ignores any x0.
    options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'Display', 'off');
    % --- x >= 0 ---
    if avoid_neg
        lb = zeros(nVoxels,1);
        ub = [];  % No upper bound.
    else
        lb = [];
        ub = [];  % No bounds.
    end
   
    % Quadratic optimization.
    [x, ~, exitflag] = quadprog(N, f, [], [], [], [], lb, ub, [], options);
    if exitflag ~= 1
        disp(['quadprog exitflag is ',num2str(exitflag),', check the covariance model.'])
    end
    
else
        
    if avoid_neg 
        
        % Unweighted NNLS.
        [~, ~, info] = nnls(At * A, At * b, struct('Accy',accy,'Order',[]));
        x  = nnls(At * A, At * b, struct('Accy',accy,'Order',info.Order));
        
    elseif ~avoid_neg
        
        % Unweighted least-squares solution without NNLS constraints.
        x  = lsqminnorm(At * A,At * b); % must faster compared with pinv
        
    end
end

% A : m (lengths in meters)
% b : kg/m^2
% x : kg/m^3

residuos = A*x - b; % kg/m^2
% hist(residuos, 20); [h, p] = kstest(residuos,'Alpha',0.05); [h, pvalue] = lillietest(residuos); % Residuals should follow a normal distribution.
% qqplot(residuos)
n = size(A,1); % Number of observations.
rmse = sqrt(sum(residuos.^2)/n) ; % kg/m^2

end
