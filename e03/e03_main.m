clear
close all
clc

format compact


%% Parameters
trg = 3;    % # of targets
p = 100;    % # of cells
q = 25;     % # of sensors
load("localization.mat")        % A, D, y


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
tic;
z0 = zeros(p+q, 1);
[z, num_iter] = ista_lasso(z0, y, G, p, q, tau, tau_Lambda, false);
t_elaps = toc;
x = z(1:p);
a = z(p+1:p+q);

% build matrix of positions
l = sqrt(p);        % side length
pos = zeros(l,l);
k = 1;
for i=l:-1:1
    for j=1:l
        pos(i,j) = x(k);
        k = k+1;
    end
end

% plot estimated targets positions
figure(1)
colormap("default")
imagesc(pos)
title("ISTA: rough targets positions")

% exploit knowledge about # of targets
x_sort = sort(x, 'descend');
smallest_accepted_value = x_sort(trg);
x_hat = x >= smallest_accepted_value;

% plot targets positions
display_CPS(x_hat, [], p, 2, "ISTA: targets positions");

% sensors under attack
a_hat = a >= 1e-3;

% print output
fprintf("ISTA\n\tTime: %f s\tIterations number: %i\n", t_elaps, num_iter);
fprintf("supp{x_hat}\n");
disp(find(x_hat));
fprintf("supp{a_hat}\n");
disp(find(a_hat));


%% k Nearest Neighbors
% inizialization
knn_min = Inf;
ind = [];
num_iter_knn = 0;       % expected nchoosek(100,3) = 161'700

% brute force approach
tic;
for i = 1:p-2
    for j = i+1:p-1
        for k = j+1:p
            temp = norm(D(:,i) + D(:,j) + D(:,k) - y, 2)^2;
            if temp < knn_min
                knn_min = temp;
                ind = [i j k];
            end
            num_iter_knn = num_iter_knn + 1;
        end
    end
end
t_elaps_knn = toc;

% create x_hat
x_hat_knn = zeros(p,1);
for i=1:trg
    x_hat_knn(ind(i)) = true;
end
x_hat_knn = x_hat_knn > 0;      % convert to logical array

% plot targets positions
display_CPS(x_hat_knn, [], p, 3, "k-NN: targets positions");

% print output
fprintf("\n\nk-NN\n\tTime: %f s\tIterations number: %i\n", t_elaps_knn, num_iter_knn);
fprintf("supp{x_hat_knn}\n");
disp(find(x_hat_knn));
fprintf("\n");