clear
close all
clc

format compact


% random noisy measurement parameters
%q_array = 10:2:29;  % [10 12 14 ... 28] --> 10 increasing values for q
q = 10;
p = 20;
k = 2;

% ista parameters
epsilon = 1e-8;
delta = 1e-12;

% simulations parameters
N_ROUTINE = 15;         % # of simulation routines
N_SIM = 20;             % # of simulations in each routine

%tau_array = linspace(0.005, 0.033, N_ROUTINE);     % to simulate with variable tau

tau = 0.015;                                        % to simulate with constant tau 
lambda_array = ...                                  %  and variable lambda
    linspace(1/(1000*tau), 1/(10*tau), N_ROUTINE);   %

% final result variables
sup_rec_cnt_array = zeros(1, N_ROUTINE);
mean_conv_time_array = zeros(1, N_ROUTINE);
min_conv_time_array = zeros(1, N_ROUTINE);
max_conv_time_array = zeros(1, N_ROUTINE);

for j = 1:N_ROUTINE
    % simulations variables
    sup_rec_cnt = 0;        % support recovery count
    num_iter_array = zeros(N_SIM, 1);
    %q = q_array(j);         % increasing q value

    % perform simulations
    for i=1:N_SIM
        [y, C, x_hat, eta] = e01_rand_noisy_mes_gen(q, p, k);
    
        % ista
        %tau = norm(C,2)^(-2) - epsilon;
        %tau = tau_array(j);
        %lambda = 1 / (100*tau);
        lambda = lambda_array(j);
        tau_lambda = tau*lambda * ones(p,1);
        z0 = zeros(p, 1);
        [x, num_iter] = ista_lasso(z0, y, C, p, 0, tau, tau_lambda, delta, false);
        
        % update vars
        if nnz(x) == nnz(x_hat) && all(find(x_hat) == find(x))
            sup_rec_cnt = sup_rec_cnt + 1;
        end
        num_iter_array(i) = num_iter;
    end
    
    sup_rec_cnt_array(j) = 100 * sup_rec_cnt/N_SIM;
    mean_conv_time_array(j) = mean(num_iter_array);
    min_conv_time_array(j) = min(num_iter_array);
    max_conv_time_array(j) = max(num_iter_array);
    
end

% configure display and plot option according to the simulation
e01_display_results;