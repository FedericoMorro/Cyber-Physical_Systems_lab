clear
close all
clc

format compact


%% Parameters
trg = 3;    % # of targets
p = 100;    % # of cells
q = 25;     % # of sensors
load("tracking_moving_targets.mat")        % A, D, Y


%% Ista
% ista parameters
lambda_1 = 10;
lambda_2 = 20;
epsilon = 1e-8;

% ista variables
G = normalize([D eye(q)]);
tau = norm(G,2)^(-2) - epsilon;
tau_Lambda = tau * [lambda_1*ones(p,1); lambda_2*ones(q,1)];

% ista

z_hat = zeros(p+q, 1);
k_fin = size(Y,2);
for k=1:k_fin
    % z^+(k) = S_tl [z_hat(k) + tau*G' * (y(k) - G*z_hat(k))]
    [z_ista, num_iter] = ista_lasso(z_hat, Y(:,k), G, p, q, tau, tau_Lambda, true);
    
    x = A * z_ista(1:p);        % x_hat(k+1) = A * x_hat^+(k)
    a = z_ista(p+1:p+q);        % a_hat(k+1) = a_hat^+(k)
    z_hat = [x; a];

    % exploit knowledge about # of targets
    x_sort = sort(x, 'descend');
    smallest_accepted_value = x_sort(trg);
    x_hat = x >= smallest_accepted_value;
    
    % plot targets positions
    display_CPS(x_hat, p, 2, "ISTA: targets positions");
    
    % sensors under attack
    a_hat = a >= 1e-1;
    
    % print output
    fprintf("supp{x_hat}\n");
    disp(find(x_hat));
    fprintf("supp{a_hat}\n");
    disp(find(a_hat));
    pause(1);
end
