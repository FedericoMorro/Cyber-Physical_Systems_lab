function [x, a, num_iter] = e02_ista_part_lasso(y, C, n, q, tau, tau_Lambda)

% parameters
delta = 1e-12;

% initialization
x = zeros(n,1);         % x(0) = 0, will store x(i)
a = zeros(q,1);
num_iter = 0;

% partial lasso variables
G = [C eye(q)];
z = [x; a];
z_next = zeros(n+q,1);    % tmp variable for x(i+1)

% iterations
exit_cond = false;
while ~exit_cond

    % gradient step
    z_grad = z + tau*G' * (y - G*z);

    % shrinkage-thresholding step
    for i=1:n+q         % element-wise
        if z_grad(i) > tau_Lambda(i)
            z_next(i) = z_grad(i) - tau_Lambda(i);
        elseif z_grad(i) < -tau_Lambda(i)
            z_next(i) = z_grad(i) + tau_Lambda(i);
        else
            z_next(i) = 0;
        end
    end

    % exit condition
    if norm(z_next-z, 2) < delta
        exit_cond = true;
    end

    % x update
    z = z_next;

    num_iter = num_iter + 1;
end

% estimated vectors
x = z(1:n);
a = z(n+1:n+q);

end