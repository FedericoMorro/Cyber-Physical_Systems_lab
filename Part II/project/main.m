clear
close all
clc

format compact


% silent: 0 plot and print, 1 print, 2 none
silent = 2;


%% Setup
% network parameters
N = 6;

Adj = zeros(N,N);
Adj(3,1) = 2;
Adj(5,1) = 2;
Adj(4,2) = 2;
Adj(6,2) = 2;

g = zeros(N,1);
g(1) = 2;
g(2) = 2;

aug_graph = [
    zeros(1,N+1)
    g       Adj
];
digr = digraph(aug_graph', {'0','1','2','3','4','5','6'});

if 1
    figure, plot(digr, 'EdgeLabel', digr.Edges.Weight)
end

% desired eigenvalues
A_des_eig = 0.5*[0 -1];         % constant
%A_des_eig = [0 0];              % ramp          
%A_des_eig = 2*[0+1i 0-1i];      % sinusoidal

% reference to be tracked
x0_ref = [1 0];

% noise
noise_vec = zeros(N,1);

% local or global observers
local_obs = false;



%% Simulation(s)
[x0_sim, y0_sim, xi_sim, yi_sim, t_sim] = ...
    coop_reg(Adj, g, A_des_eig, x0_ref, local_obs, noise_vec, silent);



%% Output elaboration
% followers step response
rise_set_time = zeros(N,2);
for i = 1:N
    s = stepinfo(yi_sim{i}, t_sim);
    rise_set_time(i,:) = [s.RiseTime, s.SettlingTime];
end
figure
plot(1:N, rise_set_time, 'o'), grid on
legend('Rise time', 'Settling time')
xlabel('Follower #'), ylabel('seconds')

% step responses plot
figure, grid on, hold on
plot(t_sim,y0_sim,'k--', 'LineWidth',1.2, 'DisplayName','S_0')
for i = 1:N
    plot(t_sim,yi_sim{i}, 'DisplayName',sprintf("S_%i", i))
end
title('Followers step response'), legend