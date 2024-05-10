function [x_mse, first_conv, a_mse] = rand_atk_2(n0, k_switch, k_fin, thr)

%n0 attacks changing sensors every k_switch iterations

% parameters
trg = 3;    % # of targets
p = 100;    % # of cells
q = 25;     % # of sensors
load("tracking_moving_targets.mat")        % A, D, Y
a_corr = 0;
a_avg = 0;
x_corr = 0;
x_mse = zeros(k_fin, 1);
a_mse = zeros(k_fin, 1);
first_conv = 0;

atk_mem = zeros(q, ceil(k_fin/k_switch));
atk_num = 1;

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

%% Time-varying attacks
% create new measurements from scratch
x_tv = zeros(p, 1);         % actual positions
pos_in = randperm(100, 3);  % random initial positions
for j = pos_in
    x_tv(j) = 1;
end

Y_tv = zeros(q, k_fin);     % measurements

atk_sensors = randperm(q, n0);        
atk_true = zeros(q, 1);         
for i=atk_sensors
    atk_true(i) = 1;
end

for i = 1:k_fin
    % compute the measurement
    Yi = D*x_tv;    

    if mod(i, k_switch) == 0
        atk_mem(:,atk_num) = atk_true;
        atk_num = atk_num + 1;
        atk_sensors = randperm(q, n0);
        atk_true = zeros(q, 1);
        for j=atk_sensors
            atk_true(j) = 1;
        end
    end

    % apply the attacks
    for j = atk_sensors(1:n0)
        Yi(j) = Yi(j)*1.5;
    end

    % update the position
    Y_tv(:,i) = Yi;
    x_tv = A*x_tv;

    
end
atk_num = 1;
atk_true = atk_mem(:,1);

% repeat the algorithm
% compute initial state to compare with estimations
z0 = zeros(p+q,1);
[z, num_iter] = ista_lasso(z0, Y_tv(:,1), G, p, q, tau, tau_lambda, delta, false);
x = z(1:p);
a = z(p+1:p+q);

x_sort = sort(x, 'descend');
smallest_accepted_value = x_sort(trg);
x = x >= smallest_accepted_value;

a_thr = thr;
a = a >= a_thr;

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
    
    % sensors under attack and error
    a_hat = a_est >= a_thr;

    if(mod(k, k_switch) == 0 && atk_num < k_fin/k_switch)
        atk_num = atk_num + 1;
        atk_true = atk_mem(:,atk_num);
    end

    a_mse(k) = norm(atk_true - a_hat, 1);

    % update actual position
    x = A*x;
    x = x > 0;      % convert to logical
    x_err = norm(x - x_hat, 1);
    if x_err == 0 && ~first_conv
        first_conv = k;
    end
    x_mse(k) = x_err;
end
