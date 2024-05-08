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


%% Preliminary analysis
% stack setups in a vector
Q_vec = {Q_4 Q_8 Q_12 Q_18};
Q_names = {'Q_4' 'Q_8' 'Q_{12}' 'Q_{18}'};

f = figure(1);
f.Position([3 4]) = [525, 400];
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
str = sprintf("ISTA");
display_CPS([], x_corr, D, [], a_corr, p, q, 2, str);
 
% print output
fprintf("\nISTA\n");
fprintf("# of iterations: %i\n", num_iter);
fprintf("supp{x}\n");
disp(supp_x_corr');
fprintf("supp{a}\n");
disp(supp_a_corr');
fprintf("\n\n");



%% DISTA

% iteration variables
num_iter = zeros(length(Q_vec),1);
k_stop = zeros(length(Q_vec),1);
x_norm_error = zeros(length(Q_vec), 20000);
a_norm_error = zeros(length(Q_vec), 20000);
k_x_conv = zeros(length(Q_vec), 1);
k_a_conv = zeros(length(Q_vec), 1);
k_x_cons = zeros(length(Q_vec), 1);
k_a_cons = zeros(length(Q_vec), 1);

% for each setup
for Q_index = 1:length(Q_vec)
    % modification of ista_lasso.m

    fprintf("DISTA - ");
    disp(cell2mat(Q_names(Q_index)))
    fprintf('Current time k: %8i', 0);

    % get current Q
    Q = cell2mat(Q_vec(Q_index));

    % initialization
    k = 1;
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

            num_iter(Q_index) = num_iter(Q_index) + 1;
        end
        
        % check if reached exit condition
        sum = 0;
        for i = 1:q
            sum = sum + norm(z_next(:,i) - z(:,i), 2)^2;
        end
        if sum < delta 
            exit_cond = true;
        end

        % variables update
        z = z_next;
        
        % check if convergence achieved
        %   by getting estimated vectors for each sensor
        x_hat = z(1:p,:);
        a_hat = z(p+1:p+q,:);

        for i = 1:q
        
            % exploit knowledge about # of targets
            x_sort = sort(x_hat(:,i), 'descend');
            smallest_accepted_value = x_sort(trg);
            x_hat(:,i) = x_hat(:,i) >= smallest_accepted_value;
            supp_x = find(x_hat(:,i));
            
            % sensors under attack
            a_hat(:,i) = a_hat(:,i) >= 0.002;
            supp_a = find(a_hat(:,i));

            % update vectors of errors
            x_norm_error(Q_index, k) = x_norm_error(Q_index, k) + norm(x_hat(:,i) - x_corr, 1);
            a_norm_error(Q_index, k) = x_norm_error(Q_index, k) + norm(a_hat(:,i) - a_corr, 1);
        end
        
        % perform average of error of each 
        x_norm_error(Q_index, k) = x_norm_error(Q_index, k) / q;
        a_norm_error(Q_index, k) = a_norm_error(Q_index, k) / q;

        % check if consesus reached
        cons_x = true;
        for i = 2:q
            if ~isequal(x_hat(:,1), x_hat(:,i))
                cons_x = false;
            end
        end
        if k > 1 && k_x_cons(Q_index) == 0 && cons_x
            k_x_cons(Q_index) = k;
        end

        cons_a = true;
        if isequal(a_hat(:,1), zeros(q,1))
            cons_a = false;
        end
        for i = 2:q
            if ~isequal(a_hat(:,1), a_hat(:,i))
                cons_a = false;
            end
        end
        if k > 1 && k_a_cons(Q_index) == 0 && cons_a
            k_a_cons(Q_index) = k;
        end

        % print
        if rem(k,100) == 0
            fprintf('\b\b\b\b\b\b\b\b%8i', k);
        end

        % update k
        k = k + 1;
    
    end

    % save k
    k_stop(Q_index) = k;

    % find convergence times
    k_x_conv(Q_index) = find(x_norm_error(Q_index,:) == 0, 1, "first");
    k_a_conv(Q_index) = find(a_norm_error(Q_index,:) == 0, 1, "first");

    % print output
    fprintf("\nTeermination condition reached in k=%i time steps, with %i total iterations, delta=%.8f\n", ...
        k, num_iter(Q_index), sum)
    fprintf("Consensus reached in: for x: %i iterations, for a: %i iterations\n", ...
        k_x_cons(Q_index), k_a_cons(Q_index))
    fprintf("Convergence to exact solution achieved in: for x: %i iterations, for a: %i iterations\n\n", ...
        k_x_conv(Q_index), k_a_conv(Q_index))
end

