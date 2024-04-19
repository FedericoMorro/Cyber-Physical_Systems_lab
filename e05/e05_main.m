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

% algorithms
G = [D eye(q)];
tau = 4e-7;
lambda = [10*ones(p,1); 0.1*ones(q,1)];
tau_lambda = tau*lambda;
delta = 1e-8;


%% Preliminary analysis
% eigenvalue analysis
eig_Q4 = maxk(abs(eig(Q_4)), 2)
eig_Q8 = maxk(abs(eig(Q_8)), 2)
eig_Q12 = maxk(abs(eig(Q_12)), 2)
eig_Q18 = maxk(abs(eig(Q_18)), 2)

% graphs plot
figure
subplot(2,2,1); plot(digraph(Q_4)); title('Q_4');
subplot(2,2,2); plot(digraph(Q_8)); title('Q_8');
subplot(2,2,3); plot(digraph(Q_12)); title('Q_{12}');
subplot(2,2,4); plot(digraph(Q_18)); title('Q_{18}');

% stack setups in a vector
Q_vec = {Q_4 Q_8 Q_12 Q_18};


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

% sensors under attack
a_sort = sort(a_est, 'descend');
smallest_accepted_value = a_sort(atk);
a_corr = a_est >= smallest_accepted_value;

% plot targets positions
str = sprintf("Correct positions");
display_CPS(x_corr, [], p, 2, str);
 
%print output
fprintf("\nISTA\n");
fprintf("# of iterations: %i\n", num_iter);
fprintf("supp{x}\n");
disp(find(x_corr));
fprintf("supp{a}\n");
disp(find(a_corr));



%% DISTA
% for each setup
for Q_cell = Q_vec
    % modification of ista_lasso.m
    Q = cell2mat(Q_cell);

    % initialization
    k = 1;
    num_iter = 0;
    z = zeros(p+q, q);      % z collects each z^(i) in R^(p+q), for each sensor
    z_next = zeros(p+q, q);     % tmp variable for z(i+1)
    
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
                elseif z_grad(s) < -tau_lambda(s)
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
    end


    %% Estimates
    % get estimated vectors for each sensor
    x_est = z(1:p);        % x_hat(k+1) = A * x_hat^+(k)
    a_est = z_ista(p+1:p+q);        % a_hat(k+1) = a_hat^+(k)
    z_hat = [x_est; a_est];

    % exploit knowledge about # of targets
    x_sort = sort(x_est, 'descend');
    smallest_accepted_value = x_sort(trg);
    x_hat = x_est >= smallest_accepted_value;
    
    % sensors under attack
    a_hat = a_est >= 1e-1;

    % update actual position
    x = A*x;
    x = x > 0;      % convert to logical

    % plot targets positions
    str = sprintf("Iteration %d", k);
    display_CPS(x_hat, x, p, 2, str);
     
    %print output
    fprintf("supp{x_hat}\n");
    disp(find(x_hat));
    fprintf("supp{a_hat}\n");
    disp(find(a_hat));

    fprintf("\nConsensus reached in %i time, with %i iterations\n", k, num_iter)
    fprintf("supp{x_hat}\n")
    fprintf("supp{a_hat}\n")

end

