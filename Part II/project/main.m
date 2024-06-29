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
g = zeros(N,1);

% binary tree configuration
% Adj(3,1) = 2;
% Adj(5,1) = 2;
% Adj(4,2) = 2;
% Adj(6,2) = 2;
% g(1) = 2;
% g(2) = 2;

% loop configuration
% Adj(3,1) = 2;
% Adj(4,3) = 2;
% Adj(2,4) = 2;
% Adj(5,2) = 2;
% Adj(6,5) = 2;
% Adj(1,6) = 2;
% g(1) = 20;
% g(2) = 20;

% linear configuration
Adj(2,1) = 1;
Adj(3,2) = 1;
Adj(4,3) = 1;
Adj(5,4) = 1;
Adj(6,5) = 1;
g(1) = 1;

aug_graph = [
    zeros(1,N+1)
    g       Adj
];
digr = digraph(aug_graph', {'0','1','2','3','4','5','6'});

if 1
    f = figure();
    f.Position([3 4]) = [525, 400];
    plot(digr, 'EdgeLabel', digr.Edges.Weight)
    title('Communication Network Topolody')
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
rise_set_time_agents = zeros(N,2);
for i = 1:N
    s = stepinfo(yi_sim{i}, t_sim, x0_ref(1)*708.27, 0, ...
        'SettlingTimeThreshold', 0.02, 'RiseTimeLimits', [0.1 0.9]);
    rise_set_time_agents(i,:) = [s.RiseTime, s.SettlingTime];
end
rise_set_time_agents
mean(rise_set_time_agents)

% rmse of response w.r.t. leader
rms_agents = zeros(N,1);
for i = 1:N
    rms_agents(i) = rms(y0_sim - yi_sim{i});
end
rms_agents
mean(rms_agents)

% plots
figure
subplot(2,1,1), plot(1:N, rise_set_time_agents, 'o'), grid on
xlim([0.5 N+0.5]), ylim([0 max(max(rise_set_time_agents))+0.5])
legend('Rise time', 'Settling time')
xlabel('Follower #'), ylabel('seconds')
title('Followers Rise and Settling Time')

subplot(2,1,2), plot(1:N, rms_agents, 'o'), grid on
xlim([0.5 N+0.5])
legend('Output RMS')
xlabel('Follower #'), ylabel('RMS')
title('Followers Output RMS')


% step responses plot
f = figure();
f.Position([3 4]) = [525, 400];
grid on, hold on
plot(t_sim,y0_sim,'k--', 'LineWidth',1, 'DisplayName','S_0')
for i = 1:N
    plot(t_sim,yi_sim{i}, 'DisplayName',sprintf("S_%i", i))
end
title('Agents step response'), legend

