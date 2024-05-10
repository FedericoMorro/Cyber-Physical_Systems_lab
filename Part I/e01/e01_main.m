clear
close all
clc

format compact

% Set simulation mode, among:
%   SimulationMode.routine_n        -> in function of routine #
%   SimulationMode.variable_q       -> in function of increasing q value
%   SimulationMode.variable_lambda  -> in function of variable lambda value
%   SimulationMode.variable_tau     -> in function of variable tau value

sim_mode = SimulationMode.routine_n;

% random noisy measurement parameters
if sim_mode == 'variable_q'
    q_array = 10:2:29;  % [10 12 14 ... 28] --> 10 increasing values for q
else
    q = 10;             % fixed q value
end

p = 20;
k = 2;

% ista parameters
epsilon = 1e-8;
delta = 1e-12;

% simulations parameters
N_ROUTINE = 10;         % # of simulation routines
N_SIM = 20;             % # of simulations in each routine

if sim_mode == 'variable_tau'
    % to simulate with variable tau
    tau_array = linspace(0.002, 0.032, N_ROUTINE);     

elseif sim_mode == 'variable_lambda'
    % to simulate with constant tau and variable lambda
    tau = 0.015;
    lambda_array = ...
        linspace(1/(1000*tau), 1/(10*tau), N_ROUTINE);
end

% final result variables
sup_rec_cnt_array = zeros(1, N_ROUTINE);
mean_conv_time_array = zeros(1, N_ROUTINE);
min_conv_time_array = zeros(1, N_ROUTINE);
max_conv_time_array = zeros(1, N_ROUTINE);

for j = 1:N_ROUTINE
    % simulations variables
    sup_rec_cnt = 0;        % support recovery count
    num_iter_array = zeros(N_SIM, 1);
    
    if sim_mode == 'variable_q'
        q = q_array(j);         % increasing q value
    end

    % perform simulations
    for i=1:N_SIM
        [y, C, x_hat, eta] = e01_rand_noisy_mes_gen(q, p, k);
        
        % ista
        if sim_mode == 'variable_tau'
            tau = tau_array(j);
            lambda = 1 / (100*tau);     % constant tau_lambda
        
        elseif sim_mode == 'variable_lambda'
            % tau defined before starting (constant)
            lambda = lambda_array(j);
        
        else 
            tau = norm(C,2)^(-2) - epsilon;
            lambda = 1 / (100*tau);
        end

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

%% display results according to chosen simulation mode

e01_display_results;