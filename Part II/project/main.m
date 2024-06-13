clear
close all
clc

format compact



%% Agents LTI parameters
A = [
    0       1
    880.87  0
];
B = [0 -9.9453]';
C = [708.27 0];
D = 0;

n = 2;
m = 1;
p = 1;

% stability analysis
eig_A = eig(A)


%% CPS paramaters
N = 6;

Adj = zeros(N,N);
Adj(2,1) = 2;
Adj(3,2) = 6;
Adj(4,3) = 1;
Adj(5,4) = 1;
Adj(6,5) = 3;

g = zeros(N,1);
g(1) = 1;

% initial conditions
x0_0 = [10 0];
xi_0 = zeros(n,1);



%% Leader controller
% static state feedback
%   u = -K*x + N*r
%   x' = (A - B*K)*x + B*N*r

% define static-state feedback loop gain
K_contr = place(A,B, 0.5*[0 -1]);           % step
%K_contr = acker(A,B, 0.5*[0 0]);            % ramp
%K_contr = place(A,B, 0.5*[0+1i 0-1i]);      % sinusoidal
A_contr = A - B*K_contr;

% subsitute A and B with controller versions
%   save olds for leader that has explicit feedback loop
A0 = A;
A = A_contr;

% impulse response (modal analysis)
impulse(ss(A,B,C,D)), grid on



%% Local controllers
% solve ARE: A'P + PA + Q - PB inv(R) B'P = 0
R = ones(m,m);
Q = eye(n);
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
out = sim('coop_reg_SVFB.slx');



%% Plot results
figure, plot(out.tout, out.x0), grid on
title('Leader S_0')
legend('x','v')

% reshape node states
x_hist = zeros(n,length(out.tout),N);
foll_n = 0;
for k = 1:n*N
    if rem(k,n) == 1, foll_n = foll_n + 1; end

    if rem(k,n) == 0
        x_hist(n,:,foll_n) = out.xi_all(k,:);
    else
        x_hist(rem(k,n),:,foll_n) = out.xi_all(k,:);
    end
end

% plot nodes behavior (only positions)
for foll_n = 1:N

    figure

    subplot(2,1,1)
    plot(out.tout, x_hist(1,:,foll_n)), hold on
    plot(out.tout, out.x0(:,1), '--'), grid on
    str = sprintf('x_%i', foll_n); legend(str,'x_0')
    str = sprintf('S_%i - position', foll_n); title(str)

    subplot(2,1,2)
    plot(out.tout, x_hist(2,:,foll_n)), hold on
    plot(out.tout, out.x0(:,2), '--'), grid on
    str = sprintf('v_%i', foll_n); legend(str,'v_0')
    str = sprintf('S_%i - velocity', foll_n); title(str)

end

% compute and plot global disagreement error
delta = zeros(n*N,length(out.tout));
delta_MS = zeros(length(out.tout));
for t = 1:length(out.tout)
    x0_bar = [out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'; out.x0(t,:)'];
    delta(:,t) = out.xi_all(:,t) - x0_bar;
    delta_MS(t) = 1/N * norm(delta(:,t))^2;
end

figure
plot(out.tout, delta_MS), grid on
title('Mean Square of global disagreement error')