clear
close all
clc

format compact


% parameters
trg = 3;    % # of targets
p = 100;    % # of cells
q = 25;     % # of sensors
load("localization.mat")        % A, D, y

% ista parameters
lambda_1 = 10;
lambda_2 = 20;
epsilon = 1e-8;

% ista variables
G = normalize([D eye(q)]);
tau = norm(G,2)^(-2) - epsilon;
tau_Lambda = tau * [lambda_1*ones(p,1); lambda_2*ones(q,1)];

% ista
[z, num_iter] = ista_lasso(y, G, p, q, tau, tau_Lambda);
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

% exploit knowledge about # of targets
x_sort = sort(x, 'descend');
smallest_accepted_value = x_sort(trg);
pos = pos >= smallest_accepted_value;

% plot targets positions
figure
colormap("default")
imagesc(pos)

% sensors under attack
a_tol = a >= 1e-3;
under_attack = find(a_tol)

% k-NN
knn = 1e10;
ind = [];

% brute force approach
for i = 1:p-2
    for j = i+1:p-1
        for k = j+1:p
            temp = norm(D(:,i) + D(:,j) + D(:,k) - y)^2;
            if temp < knn
                knn = temp;
                ind = [i j k];
            end
        end
    end
end

% translate vector indeces to matrix' and plot
pos_knn = zeros(l, l);
for i = 1:3
    pos_knn(l - floor(ind(i)/l), rem(ind(i), l)) = 1;
end
figure
colormap("default")
imagesc(pos_knn)
