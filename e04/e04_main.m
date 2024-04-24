clear
close all
clc

format compact


% parameters
trg = 3;    % # of targets
atk = 2;    % # of attacks
p = 100;    % # of cells
q = 25;     % # of sensors
load("tracking_moving_targets.mat")        % A, D, Y


%% Ista
% ista parameters
lambda_1 = 10;
lambda_2 = 20;
epsilon = 1e-8;
delta = 1e-12;

% ista variables
G = normalize([D eye(q)]);
tau = norm(G,2)^(-2) - epsilon;
tau_lambda = tau * [lambda_1*ones(p,1); lambda_2*ones(q,1)];

% compute initial state to compare with estimations
z0 = zeros(p+q,1);
[z, num_iter] = ista_lasso(z0, Y(:,1), G, p, q, tau, tau_lambda, delta, false);
x = z(1:p);
a = z(p+1:p+q);

x_sort = sort(x, 'descend');
smallest_accepted_value = x_sort(trg);
x = x >= smallest_accepted_value;

a_sort = sort(a, 'descend');
smallest_accepted_value = a_sort(atk);
a = a >= smallest_accepted_value;

% dynamic sse ista
z_hat = zeros(p+q,1);
k_fin = size(Y,2);
for k=1:k_fin
    % z^+(k) = S_tl [z_hat(k) + tau*G' * (y(k) - G*z_hat(k))]
    [z_ista, ~] = ista_lasso(z_hat, Y(:,k), G, p, q, tau, tau_lambda, delta, true);
    
    x_est = A * z_ista(1:p);        % x_hat(k+1) = A * x_hat^+(k)
    a_est = z_ista(p+1:p+q);        % a_hat(k+1) = a_hat^+(k)
    z_hat = [x_est; a_est];

    % exploit knowledge about # of targets
    x_sort = sort(x_est, 'descend');
    smallest_accepted_value = x_sort(trg);
    x_hat = x_est >= smallest_accepted_value;
    
    % sensors under attack
    a_hat = a_est >= 1;

    % update actual position
    x = A*x;
    x = x > 0;      % convert to logical

    % plot targets positions
    str = sprintf("Iteration %d", k);
    display_CPS(x_hat, x, D, a_hat, a, p, q, 2, str);

    pause(0.5)
end




%% Time-varying attacks
% create new measurements from scratch
x_tv = zeros(p, 1);         % actual positions
pos_in = [13, 40, 84];   
for j = pos_in
    x_tv(j) = 1;
end

Y_tv = zeros(q, k_fin);     % measurements

n_atks = 4;
atk_sensors = randperm(q, n_atks);     % sensors under attack
n0 = 2;                     % initial number of attacked sensors
k_switch = 24;              % time after which the second attack starts

for i = 1:k_fin
    % compute the measurement
    Yi = D*x_tv;
    
    % apply the attacks
    for j = atk_sensors(1:n0)
        Yi(j) = Yi(j)*1.5;
    end

    if i > k_switch         % second wave of attacks
        for j = atk_sensors(n0+1:n_atks)
            Yi(j) = Yi(j)*1.5;
        end
    end

    % update the position
    Y_tv(:,i) = Yi;
    x_tv = A*x_tv;
end

% repeat the algorithm
% compute initial state to compare with estimations
z0 = zeros(p+q,1);
[z, num_iter] = ista_lasso(z0, Y_tv(:,1), G, p, q, tau, tau_lambda, delta, false);
x = z(1:p);
a = z(p+1:p+q);

x_sort = sort(x, 'descend');
smallest_accepted_value = x_sort(trg);
x = x >= smallest_accepted_value;

a = a >= 1e-1;

% dynamic sse ista
z_hat = zeros(p+q,1);
k_fin = size(Y_tv,2);
for k=1:k_fin
    % z^+(k) = S_tl [z_hat(k) + tau*G' * (y(k) - G*z_hat(k))]
    [z_ista, ~] = ista_lasso(z_hat, Y_tv(:,k), G, p, q, tau, tau_lambda, delta, true);
    
    x_est = A * z_ista(1:p);        % x_hat(k+1) = A * x_hat^+(k)
    a_est = z_ista(p+1:p+q);        % a_hat(k+1) = a_hat^+(k)
    z_hat = [x_est; a_est];

    % exploit knowledge about # of targets
    x_sort = sort(x_est, 'descend');
    smallest_accepted_value = x_sort(trg);
    x_hat = x_est >= smallest_accepted_value;
    
    % sensors under attack
    a_hat = a_est >= 1;

    % update actual position
    x = A*x;
    x = x > 0;      % convert to logical

    % plot targets positions
    str = sprintf("Iteration %d", k);
    display_CPS(x_hat, x, D, a_hat, a, p, q, 3, str);

    pause(0.5);
end