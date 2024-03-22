function [x, num_iter] = e01_ista_lasso(y, C, p, q, tau, lambda)

% parameters
gamma = tau*lambda;
delta = 1e-12;

% initialization
x = zeros(p,1);         % x(0) = 0, will store x(i)
x_next = zeros(p,1);    % tmp variable for x(i+1)
num_iter = 0;

% iterations
exit_cond = false;
while ~exit_cond

    % gradient step
    x_grad = x + tau*C' * (y - C*x);

    % shrinkage-thresholding step
    for i=1:p       % element-wise
        if x_grad(i) > gamma
            x_next(i) = x_grad(i) - gamma;
        elseif x_grad(i) < -gamma
            x_next(i) = x_grad(i) + gamma;
        else
            x_next(i) = 0;
        end
    end

    % exit condition
    if norm(x_next-x, 2) < delta
        exit_cond = true;
    end

    % x update
    x = x_next;

    num_iter = num_iter + 1;
end

end