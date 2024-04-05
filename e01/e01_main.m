clear
close all
clc

format compact


% random noisy measurement parameters
q = 10;
p = 20;
k = 2;

% ista parameters
epsilon = 1e-8;

% simulations parameters
N_SIM = 20;

% simulations variables
sup_rec_cnt = 0;        % support recovery count
num_iter_array = zeros(N_SIM, 1);

% perform simulations
for i=1:N_SIM
    [y, C, x_hat, eta] = e01_rand_noisy_mes_gen(q, p, k);

    % ista
    tau = norm(C,2)^(-2) - epsilon;
    lambda = 1 / (100*tau);
    tau_Lambda = tau*lambda * ones(p,1);
    [x, num_iter] = ista_lasso(y, C, p, 0, tau, tau_Lambda);
    
    % update vars
    if nnz(x) == nnz(x_hat)
        sup_rec_cnt = sup_rec_cnt + 1;
    end
    num_iter_array(i) = num_iter;
end

% print results
fprintf("Support recovery rate\n\t%i%%\n", sup_rec_cnt/N_SIM * 100);
fprintf("Convergence time analysis\n");
fprintf("\tMean: %.2f\n\tMin: %i\n\tMax: %i\n", ...
    mean(num_iter_array), min(num_iter_array), max(num_iter_array));