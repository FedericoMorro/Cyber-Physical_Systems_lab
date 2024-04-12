function [z, num_iter] = ista_lasso(z0, y, G, p, q, tau, tau_Lambda, one_iter)

% parameters
delta = 1e-12;

% initialization
num_iter = 0;
z = z0;
z_next = zeros(p+q, 1);     % tmp variable for z(i+1)

% iterations
exit_cond = false;
while ~exit_cond

    % gradient step
    z_grad = z + tau*G' * (y - G*z);

    % shrinkage-thresholding step
    for i=1:p+q         % element-wise
        if z_grad(i) > tau_Lambda(i)
            z_next(i) = z_grad(i) - tau_Lambda(i);
        elseif z_grad(i) < -tau_Lambda(i)
            z_next(i) = z_grad(i) + tau_Lambda(i);
        else
            z_next(i) = 0;
        end
    end

    % exit condition
    if norm(z_next-z, 2) < delta || one_iter
        exit_cond = true;
    end

    % x update
    z = z_next;

    num_iter = num_iter + 1;
end

end