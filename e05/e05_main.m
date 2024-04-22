clear
close all
clc

format compact


%% Parameters
% CPS
p = 100;    % # of cells
q = 25;     % # of sensors
trg = 2;    % # of targets
atk = 2;    % # of sensors under attack
load('distributed_localization_data.mat');

% DISTA algorithm
G = [D eye(q)];
tau = 4e-7;
lambda = [10*ones(p,1); 0.1*ones(q,1)];
tau_lambda = tau*lambda;
delta = 1e-8;
MAX_K = 2e5;


%% Preliminary analysis
% stack setups in a vector
Q_vec = {Q_4 Q_8 Q_12 Q_18};
Q_names = {'Q_4' 'Q_8' 'Q_{12}' 'Q_{18}'};

figure(1)
for Q_index = 1:length(Q_vec)
    % get current Q
    Q = cell2mat(Q_vec(Q_index));

    % effective spectral radius analysis
    fprintf("ESR of ");
    disp(cell2mat(Q_names(Q_index)));
    disp(maxk(abs(eig(Q)),2)');

    % plot graph in subplot
    subplot(2,2,Q_index);
    plot(digraph(Q));
    title(Q_names(Q_index));
end


%% ISTA
% get correct solution with ISTA
z0 = zeros(p+q,1);
G_ista = normalize(G);
tau_ista = norm(G_ista,2)^(-2) - 1e-8;
tau_lambda_ista = tau * [10*ones(p,1); 0.1*ones(q,1)];
delta_ista = 1e-5;
[z_corr, num_iter] = ista_lasso(z0, y, G_ista, p, q, tau_ista, tau_lambda_ista, delta_ista, false);

% estimates
x_est = z_corr(1:p);
a_est = z_corr(p+1:p+q);

% exploit knowledge about # of targets
x_sort = sort(x_est, 'descend');
smallest_accepted_value = x_sort(trg);
x_corr = x_est >= smallest_accepted_value;
supp_x_corr = find(x_corr);

% sensors under attack
a_sort = sort(a_est, 'descend');
smallest_accepted_value = a_sort(atk);
a_corr = a_est >= smallest_accepted_value;
supp_a_corr = find(a_corr);

% plot targets positions
str = sprintf("Correct positions");
display_CPS(x_corr, [], p, 2, str);
 
%print output
fprintf("\nISTA\n");
fprintf("# of iterations: %i\n", num_iter);
fprintf("supp{x}\n");
disp(supp_x_corr');
fprintf("supp{a}\n");
disp(supp_a_corr');
fprintf("\n");



%% DISTA
% for each setup
for Q_index = 1:length(Q_vec)
    % modification of ista_lasso.m

    fprintf("\nDISTA - ");
    disp(cell2mat(Q_names(Q_index)))
    fprintf('Current time k: %8i', 0);

    % get current Q
    Q = cell2mat(Q_vec(Q_index));

    % initialization
    k = 1;
    num_iter = 0;
    z = zeros(p+q, q);      % z collects each z^(i) in R^(p+q), for each sensor
    z_next = zeros(p+q, q);     % tmp variable for z_{k+1}
    
    % iterations
    exit_cond = false;
    while ~exit_cond
    
        % for each sensor
        for i = 1:q
    
            % compute Q contributions
            z_grad = 0;
            for j = 1:q
                z_grad = z_grad + Q(i,j) * z(:,j);
            end

            % compute gradient step
            z_grad = z_grad + tau*G(i,:)' * (y(i) - G(i,:)*z(:,i));
    
            % apply shrinkage-thresholding operator
            for s = 1:p+q         % element-wise on i-th column of z <-> i-th sensor
                if z_grad(s) > tau_lambda(s)
                    z_next(s,i) = z_grad(s) - tau_lambda(s);
                elseif z_grad(s) < - tau_lambda(s)
                    z_next(s,i) = z_grad(s) + tau_lambda(s);
                else
                    z_next(s,i) = 0;
                end
            end

            num_iter = num_iter + 1;
        end
        
        % check if reached exit condition
        sum = 0;
        for i = 1:q
            sum = sum + norm(z_next(:,i) - z(:,i));
        end
        if sum < delta 
            exit_cond = true;
        end
    
        % variables update
        z = z_next;
        k = k + 1;

        if rem(k,10000) == 0
            fprintf('\b\b\b\b\b\b\b\b%8i', k);
        end
        if k >= MAX_K
            exit_cond = true;
        end
    end


    % get estimated vectors for each sensor
    n_corr_x = 0;
    n_corr_a = 0;
    for i = 1:q
        x_est = z(1:p, i);
        a_est = z(p+1:p+q, i);
    
        % exploit knowledge about # of targets
        x_sort = sort(x_est, 'descend');
        smallest_accepted_value = x_sort(trg);
        x_hat = x_est >= smallest_accepted_value;
        supp_x = find(x_hat);
        
        % sensors under attack
        a_hat = a_est >= 0.002;
        supp_a = find(a_hat);

        % count sensors that have correct estimations
        if isequal(supp_x(1:length(supp_x_corr)), supp_x_corr)
            n_corr_x = n_corr_x + 1;
        end
        if isequal(supp_a(1:length(supp_a_corr)), supp_a_corr)
            n_corr_a = n_corr_a + 1;
        end
    
        % plot targets positions
        %display_CPS(x_hat, x_corr, p, 2, str);
    end
     
    % print output
    fprintf("\nConsensus reached in k=%i time steps, with %i total iterations, delta=%.9f\n", ...
        k, num_iter, sum)
    fprintf("Number of sensor with correct supports: %i/%i for x, %i/%i for a\n", ...
        n_corr_x, q, n_corr_a, q)
end

