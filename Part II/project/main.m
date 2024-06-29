clear all
close all
clc

format compact


% silent: 0 plot and print, 1 print, 2 none
silent = 2;

%% Setup
% network parameters
N = 6;

[Adj, g] = network_config("", N, 0);

ref_type = "step";
A_des_eig = reference_config(ref_type);

% reference to be tracked
x0_ref = [1 0];
output_fact = 708.27;

% noise
noise_vec = zeros(N+1,1);
%noise_vec(6) = 1;

% local or global observers
local_obs = false;

% algorithm parameters
n = 2; m = 1;
par.c_fact = 2;
par.co_fact = 2;
par.R = eye(m);
par.Q = eye(n);
par.Ro = eye(m);
par.Qo = eye(n);



%% Simulation(s)
[x0_sim, y0_sim, xi_sim, yi_sim, ui_sim, t_sim] = ...
    coop_reg(Adj, g, A_des_eig, x0_ref, par, local_obs, noise_vec, silent);



%% Output elaboration

% rmse of response w.r.t. leader
rms_agents = zeros(N,1);
for i = 1:N
    rms_agents(i) = rms(y0_sim - yi_sim{i});
end
rms_agents
mean(rms_agents)

% command inputs norm
effort_agents = zeros(N,1);
for i = 1:N
    effort_agents(i) = norm(ui_sim{i});
end
effort_agents
mean(effort_agents)

% global disagreement error
delta = zeros(n*N, length(x0_sim));
delta_MS = zeros(length(x0_sim),1);
for t = 1:length(x0_sim)
    x0_bar = kron(ones(N,1), x0_sim(t,:)');

    xi_all = [];
    for i = 1:N
        xi_all = [xi_all; xi_sim{i}(:,t)];
    end

    delta(:,t) = xi_all - x0_bar;
    delta_MS(t) = 1/N * norm(delta(:,t))^2;
end

% plots
% output responses plot
f = figure();
f.Position([3 4]) = [525, 400];
grid on, hold on
plot(t_sim,y0_sim,'k--', 'LineWidth',1, 'DisplayName','S_0')
for i = 1:N
    plot(t_sim,yi_sim{i}, 'DisplayName',sprintf("S_%i", i))
end
title('Agents step response'), legend

% command inputs plot
f = figure();
f.Position([3 4]) = [525, 400];
grid on, hold on
for i = 1:N
    plot(t_sim,ui_sim{i}, 'DisplayName',sprintf("S_%i", i))
end
title('Agents command inputs'), legend

% global disagreement error
figure
plot(t_sim, delta_MS), grid on
title('Mean Square of global disagreement error')


%% Time analysis
if ref_type == "step"
    % followers step response
    times_agents = zeros(N,2);
    for i = 1:N
        s = stepinfo(yi_sim{i}, t_sim, x0_ref(1)*output_fact, 0, ...
            'SettlingTimeThreshold', 0.02, 'RiseTimeLimits', [0.1 0.9]);
        times_agents(i,:) = [s.SettlingTime, s.RiseTime];
    end
else
    % followers settling time
    times_agents = zeros(N,1);
    for i = 1:N
        for t = flip(t_sim)
            if abs(yi_sim{i}(t) - y0_sim(t)) > x0_ref(1)*output_fact*0.02
                times_agents(i) = t;
                break;
            end
        end
    end
end
times_agents
mean(times_agents)

% plots
f = figure();
f.Position([1 2 3 4]) = [0, 0, 525, 2*400];

subplot(3,1,1), plot(1:N, times_agents, 'o'), grid on
xlim([0.5 N+0.5]), ylim([0 max(max(times_agents))+0.5])
if ref_type == "step"
    legend('Settling time', 'Rise time'), title('Followers Rise and Settling Time')
else
    legend('Settling time'), title('Followers Settling Time')
end
xlabel('Follower #'), ylabel('seconds')

subplot(3,1,2), plot(1:N, rms_agents, 'o'), grid on
xlim([0.5 N+0.5]), ylim([0.95*min(rms_agents) 1.05*max(rms_agents)])
legend('Output RMS')
xlabel('Follower #'), ylabel('RMS')
title('Followers Output RMS')

subplot(3,1,3), plot(1:N, effort_agents, 'o'), grid on
xlim([0.5 N+0.5]), ylim([0.95*min(effort_agents) 1.05*max(effort_agents)])
legend('Command effort norm')
xlabel('Follower #'), ylabel('norm')
title('Followers Command effort')