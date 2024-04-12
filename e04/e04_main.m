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
    [z_ista, num_iter] = ista_lasso(z_hat, Y(:,k), G, p, q, tau, tau_Lambda, true);
    
    x = A*z_ista(1:p);
    a = z_ista(p+1:p+q);
    z_hat = [x; a];
end
