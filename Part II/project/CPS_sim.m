function [metrics] = CPS_sim(ref_type, par, local_obs, silent)



%% Setup
% network parameters
N = 6;
n = 2;

[Adj, g] = network_config("", N, silent);

% reference to be tracked
[A_des_eig, x0_ref, t_fin] = reference_config(ref_type, 1, 0.5);
output_fact = 708.27;

% noise
noise_vec = zeros(N+1,1);


%% Simulation(s)
[x0_sim, y0_sim, xi_sim, yi_sim, yt_sim, ui_sim, ei_sim, t_sim] = ...
    coop_reg(Adj, g, A_des_eig, x0_ref, par, local_obs, noise_vec, t_fin, silent);



%% Output elaboration

sim_len = length(t_sim);

% rmse of response w.r.t. leader
rms_output_agents = zeros(N,1);
for i = 1:N
    rms_output_agents(i) = rms(y0_sim - yi_sim{i});
end
rms_output_agents
mean_rms_output_agents = mean(rms_output_agents);

% command inputs norm
effort_agents = zeros(N,1);
for i = 1:N
    effort_agents(i) = norm(ui_sim{i})^2;
end
effort_agents
mean_effort_agents = mean(effort_agents);

% global disagreement error
delta = zeros(n*N, sim_len);
delta_RMS = zeros(sim_len,1);
for t_ind = 1:sim_len
    x0_bar = kron(ones(N,1), x0_sim(t_ind,:)');

    xi_all = [];
    for i = 1:N
        xi_all = [xi_all; xi_sim{i}(:,t_ind)];
    end

    delta(:,t_ind) = xi_all - x0_bar;
    delta_RMS(t_ind) = rms(delta(:,t_ind));
end
mean_delta_RMS = mean(delta_RMS);

% time of zero disagreement error
t_ga = -1;
for t_ind = sim_len:-1:1
    if delta_RMS(t_ind) > 1e-2
        t_ga = t_sim(t_ind);
        break;
    end
end
t_ga

% average y_tilde
yt_avg = zeros(sim_len,1);
for i = 1:N
    yt_avg = yt_avg + abs(yt_sim{i}(:));
end
yt_avg = yt_avg / N;
mean_yt_avg = mean(yt_avg);

% observer estimation error (average of states)
avg_obs_err = {N};
for i = 1:N
    avg_obs_err{i} = zeros(sim_len,1);
    for t_ind = 1:sim_len
        avg_obs_err{i}(t_ind) = norm(ei_sim{i}(t_ind,:));
    end
end
% average
avg_avg_obs_err = zeros(sim_len, 1);
for i = 1:N
    avg_avg_obs_err = avg_avg_obs_err + avg_obs_err{i}(:);
end
avg_avg_obs_err = avg_avg_obs_err / N;
mean_avg_avg_obs_err = mean(avg_avg_obs_err);



%% Plots
if silent < 2

    % output responses plot
    f = figure();
    f.Position([3 4]) = [525, 400];
    grid on, hold on
    plot(t_sim,y0_sim,'k--', 'LineWidth',1, 'DisplayName','S_0')
    for i = 1:N
        plot(t_sim,yi_sim{i}, 'DisplayName',sprintf("S_%i", i))
    end
    title('Agents Output Response'), legend
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
    
    % global disagreement error
    figure
    plot(t_sim, delta_RMS), grid on
    title('RMS of global disagreement error')
    xlabel('Time [s]'), ylabel('MS(x_{all} - x_{bar})')
    
    % y_tilde and observer error
    figure
    subplot(2,1,1), hold on
    plot(t_sim,yt_avg, 'k--', 'LineWidth',1, 'DisplayName', 'Avg')
    for i = 1:N
        plot(t_sim,yt_sim{i}, 'DisplayName',sprintf("S_%i", i))
    end
    title('Agents tracking error'), legend(), grid on
    xlabel('Time [s]'), ylabel('y_{tilde}')
    
    subplot(2,1,2), hold on
    plot(t_sim,avg_avg_obs_err, 'k--', 'LineWidth',1, 'DisplayName', 'Avg')
    for i = 1:N
        plot(t_sim,avg_obs_err{i}, 'DisplayName',sprintf("S_%i", i))
    end
    title('Agents observer error'), legend(), grid on
    xlabel('Time [s]'), ylabel('avg(x - x_{hat})')

end



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
        for t_ind = sim_len:-1:1
            if abs(yi_sim{i}(t_ind) - y0_sim(t_ind)) > max(x0_ref)*output_fact*0.02
                times_agents(i) = t_sim(t_ind);
                break;
            end
        end
    end
end
times_agents
mean_times_agents = mean(times_agents);

if silent < 2
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
    
    subplot(3,1,2), plot(1:N, rms_output_agents, 'o'), grid on
    xlim([0.5 N+0.5]), ylim([0.95*min(rms_output_agents) 1.05*max(rms_output_agents)])
    legend('Output RMS')
    xlabel('Follower #'), ylabel('RMS')
    title('Followers Output Error RMS')
    
    subplot(3,1,3), plot(1:N, effort_agents, 'o'), grid on
    xlim([0.5 N+0.5]), ylim([0.95*min(effort_agents) 1.05*max(effort_agents)])
    legend('Command effort energy')
    xlabel('Follower #'), ylabel('Energy ||u||^2')
    title('Followers Command Effort')

end



%% Distance of agents

% farthest agent from others
dist_vec = zeros(sim_len,1);
for t_ind = 1:sim_len
    dist_matr = zeros(N,N);
    for i = 1:N
        for j = 1:N
            dist_matr(i,j) = norm(yi_sim{i}(t_ind) - yi_sim{j}(t_ind));
        end
    end
    dist_vec(t_ind) = norm(dist_matr, 1);
end
mean_dist_vec = mean(dist_vec);

% plot
if silent < 2
    figure
    plot(t_sim, dist_vec), grid on
    title('Maximum Distance of Agents Outputs')
    xlabel('Time [s]'), ylabel('dist')
end



%% Display simulation numerical results

mean_delta_RMS
mean_yt_avg

disp('--------------------')

mean_rms_output_agents
t_ga
mean_effort_agents
mean_times_agents
mean_dist_vec
mean_avg_avg_obs_err

metrics = [
    mean_rms_output_agents
    t_ga
    mean_effort_agents
    mean_times_agents(1)
    mean_dist_vec
    mean_avg_avg_obs_err
];