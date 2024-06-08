clear
close all
clc

format compact



%% Setup
% physical parameters
m1 = 1.1;
m2 = 0.9;
k1 = 1.5;
k2 = 1;

% xi = [posi_1 veli_1 posi_2 veli_2]
% xi_dot = [
%   posi_1_dot = veli_1
%   veli_1_dot = 1/m1 * (-k1*posi_1 + k2*(posi_2 - posi_1) + ui
%   posi_2_dot = veli_2
%   veli_2_dot = 1/m2 * (-k2*(posi_2 - posi_1)
% ]

% system variables
A = [
    0               1       0           0
    (-k1-k2)/m1     0       k2/m1       0
    0               0       0           1
    k2/m2           0       -k2/m2      0
];
B = [
    0
    1/m1
    0
    0
];
C = eye(4);
D = zeros(4,1);

Adj = zeros(6,6);
Adj(2,1) = 2;
Adj(3,2) = 6;
Adj(4,3) = 1;
Adj(5,4) = 1;
Adj(6,5) = 3;

g = zeros(6,1);
g(1) = 1;

% initial conditions
x0_0 = [10 0 0 0];
xi_0 = zeros(4,1);



%% Leader controller
% static state feedback
%   u = -K*x + N*r
%   x' = (A - B*K)*x + B*N*r

% define static-state feedback loop gain
K_contr = place(A,B, 0.5*[-1 -2 -3 -4]);
A_contr = A - B*K_contr;

% define static-state feedback reference gain
dc = dcgain(ss(A_contr,B,C,D));
for k = 1:length(dc)
    if dc(k) == 0; N_contr(k) = 1; else N_contr(k) = 1/dc(k); end
end
B_contr = B.*N_contr';

% reference signal
ref_const = [1 0 0 0];

% subsitute A and B with controller versions
%   save olds for leader that has explicit feedback loop
A0 = A;
B0 = B;
A = A_contr;
B = B_contr;



%% Local controllers
% solve ARE: A'P + PA + Q - PB inv(R) B'P = 0
R = 1;
Q = eye(4);
P = are(A, B*inv(R)*B', Q);

% compute K
K = inv(R)*B'*P;

% compute coupling gain c
in_deg = diag(ones(6,1));    % in-degrees
L = in_deg - Adj;            % Laplacian matrix of graph
G = diag(g);            % pinning matrix
lambda = eig(L+G);       
c_min = 1 / (2*min(real(lambda)));
c = 2*c_min;



%% Simulation
t_sim = 50;
out = sim('multi_agent_SBVF.slx');

% plot results
figure, plot(out.tout, out.x0), grid on
title('Leader S_0')
legend('x_1','v_1','x_2','v_2')

% reshape node states
x_hist = zeros(4,length(out.tout),6);
foll_n = 0;
for k = 1:24
    if rem(k,4) == 1, foll_n = foll_n + 1; end

    if rem(k,4) == 0
        x_hist(4,:,foll_n) = out.xi_all(k,:);
    else
        x_hist(rem(k,4),:,foll_n) = out.xi_all(k,:);
    end
end

% plot nodes behavior (only positions)
for foll_n = 1:6
    str = sprintf('Follower S_%i', foll_n);
    figure, hold on, title(str)

    plot(out.tout, x_hist(1,:,foll_n), 'b')
    plot(out.tout, out.x0(:,1), 'b--')
    plot(out.tout, x_hist(3,:,foll_n), 'r')
    plot(out.tout, out.x0(:,3), 'r--')

    grid on
end

% compute and plot global disagreement error
delta = zeros(24,length(out.tout));
delta_MS = zeros(length(out.tout));
for t = 1:length(out.tout)
    x0_bar = [out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'];
    delta(:,t) = out.xi_all(:,t) - x0_bar;
    delta_MS(t) = 1/6 * norm(delta(:,t))^2;
end

figure
plot(out.tout, delta_MS), grid on
title('Mean Square of global disagreement error')