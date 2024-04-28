clear
close all
clc

format compact


% random noisy measurement parameters
n = 10;
q = 20;
h = 2;
aware = false;

% ista parameters
epsilon = 1e-8;
tau_lambda_fact = 2e-3;
tau_lambda = tau_lambda_fact * [zeros(n,1); ones(q,1)];
delta = 1e-12;

% simulations parameters
N_ROUTINE = 10;
N_SIM = 20;

% final result variables
att_detection_array = zeros(1, N_ROUTINE);
max_estim_acc_array = zeros(1, N_ROUTINE);
mean_estim_acc_array = zeros(1, N_ROUTINE);
min_estim_acc_array = zeros(1, N_ROUTINE);

for j=1:N_ROUTINE
    % simulations variables
    att_det_cnt = 0;         % attack detction count
    estim_acc_array = zeros(N_SIM, 1);      % estimation accuracy
    
    % perform simulations
    for i=1:N_SIM
        [y, C, x_hat, a_hat, eta] = e02_rand_noisy_mes_gen(n, q, h, aware);
    
        % ista
        tau = norm(C,2)^(-2) - epsilon;
        G = [C eye(q)];
        z0 = zeros(n+q, 1);
        [z, num_iter] = ista_lasso(z0, y, G, n, q, tau, tau_lambda, delta, false);
        
        % estimated vectors
        x = z(1:n);
        a = z(n+1:n+q);
        
        % update variables
        if nnz(a_hat) == nnz(a) && all(find(a_hat) == find(a))
            att_det_cnt = att_det_cnt + 1;
        end
        estim_acc_array(i) = norm(x-x_hat, 2)^2;
    end
    
    % print results
    %fprintf("Attack detection rate\n\t%i%%\n", 100*att_det_cnt/N_SIM );
    %fprintf("Average estimation accuracy\n\t%.2f\n", mean(estim_acc_array));

    att_detection_array(j) = 100*att_det_cnt/N_SIM;
    max_estim_acc_array(j) = max(estim_acc_array); 
    mean_estim_acc_array(j) = mean(estim_acc_array);
    min_estim_acc_array(j) = min(estim_acc_array);
end

e02_display_results;