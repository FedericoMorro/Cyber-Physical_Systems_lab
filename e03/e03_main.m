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
[z, num_iter] = ista_lasso(y, G, p, q, tau, tau_Lambda);
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
figure
colormap("default")
imagesc(pos)
title("ISTA: rough targets positions")

% exploit knowledge about # of targets
x_sort = sort(x, 'descend');
smallest_accepted_value = x_sort(trg);
pos = pos >= smallest_accepted_value;

% plot targets positions
figure
colormap("default")
imagesc(pos)
title("ISTA: targets positions")

% sensors under attack
a_tol = a >= 1e-3;
under_attack = find(a_tol);

% print output
fprintf("ISTA\n\tTime: %f s\tIterations number: %i\n", t_elaps, num_iter);
fprintf("Targets positions\n");
disp(find(pos));
fprintf("Indexes sensors under attack\n");
disp(under_attack);


%% k Nearest Neighbors
% inizialization
knn_min = Inf;
ind = [];
num_iter_knn = 0;       % expected nchoosek(100,3) = 161'700

% brute force approach
tic
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

% translate vector indices to matrix' and plot
pos_knn = zeros(l,l);
for i = 1:trg
    pos_knn(l - floor(ind(i)/l), rem(ind(i),l)) = 1;
end
figure
colormap("default")
imagesc(pos_knn)
title("k-NN: targets positions")

% print output
fprintf("\n\nk-NN\n\tTime: %f s\tIterations number: %i\n", t_elaps_knn, num_iter_knn);
fprintf("Targets positions\n");
disp(find(pos_knn));
fprintf("\n");