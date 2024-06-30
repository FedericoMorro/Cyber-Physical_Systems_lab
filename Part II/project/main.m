clear
close all
clc

format compact


% silent: 0 plot and print, 1 print, 2 none
%   -1 for (also) impulse response
silent = 2;

%% Setup
% network parameters
N = 6;

[Adj, g] = network_config("", N, 0);

% reference to be tracked
ref_type = "sinusoidal";
[A_des_eig, x0_ref, t_fin] = reference_config(ref_type, 1, 0.5);
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
[x0_sim, y0_sim, xi_sim, yi_sim, yt_sim, ui_sim, t_sim] = ...
    coop_reg(Adj, g, A_des_eig, x0_ref, par, local_obs, noise_vec, t_fin, silent);



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

% average y_tilde
yt_avg = zeros(length(t_sim),1);
for i = 1:N
    yt_avg = yt_avg + yt_sim{i}(:);
end
yt_avg = yt_avg / N;

% plots
% output responses plot
f = figure();
f.Position([3 4]) = [525, 400];
grid on, hold on
plot(t_sim,y0_sim,'k--', 'LineWidth',1, 'DisplayName','S_0')
for i = 1:N
    plot(t_sim,yi_sim{i}, 'DisplayName',sprintf("S_%i", i))
end
title('Agents output response'), legend
xlabel('Time [s]'), ylabel('y_i')

% command inputs plot
f = figure();
f.Position([3 4]) = [525, 400];
grid on, hold on
for i = 1:N
    plot(t_sim,ui_sim{i}, 'DisplayName',sprintf("S_%i", i))
end
title('Agents command inputs'), legend
xlabel('Time [s]'), ylabel('u_i')

% global disagreement error and y_tilde
figure
subplot(2,1,1), plot(t_sim, delta_MS), grid on
title('Mean Square of global disagreement error')
xlabel('Time [s]'), ylabel('MS(x_{all} - x_{bar})')

subplot(2,1,2), hold on
plot(t_sim,yt_avg, 'k--', 'LineWidth',1, 'DisplayName', 'Avg')
for i = 1:N
    plot(t_sim,yt_sim{i}, 'DisplayName',sprintf("S_%i", i))
end
%title('Agents tracking error'), legend, grid on
xlabel('Time [s]'), ylabel('y_{tilde}')



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
        for t_ind = length(t_sim):-1:1
            if abs(yi_sim{i}(t_ind) - y0_sim(t_ind)) > max(x0_ref)*output_fact*0.02
                times_agents(i) = t_sim(t_ind);
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