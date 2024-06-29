function [x0_sim, y0_sim, xi_sim, yi_sim, ui_sim, t_sim] = ...
    coop_reg(Adj, g, A_des_eig, x0_ref, local_obs, noise_vec, silent)
% Cooperative regulation problem

% Adj: network adjacency matrix
% g: connections to leader node
% A_des_eig: desired eigenvalues of local controllers
% x0_ref: set initial state of x0 which will change regulation
% local_obs: true/false, use local or global observers
% noise_vec: vec[0/1], not noise or noise on i-th output
% silent: 0 plot and print, 1 print, 2 none


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
if silent < 2
    eig_A = eig(A)'
end



%% CPS paramaters
N = 6;

% in-degree
in_deg = zeros(1,N);
for i = 1:N
    in_deg(i) = nnz(Adj(i,:)) + nnz(g(i));
end

% 0 if no noise, 1 if noise
noise_mean = 0;     noise_var = 10;   noise_seed = randi(100);
leader_noise = 0;



%% Leader Controller
% static state feedback
%   u = -K*x
%   x' = (A - B*K)*x

% define static-state feedback loop gain
if A_des_eig(1) ~= A_des_eig(2)
    K_contr = place(A,B, A_des_eig);
else
    K_contr = acker(A,B, A_des_eig);
end
A_contr = A - B*K_contr;

if silent < 2
    eig_Acontr = eig(A_contr)'
end

% reference signal gain for x1
ss_ref = 1;

s = tf('s');
sys = minreal(zpk(inv(s*eye(n)-A_contr)));
ss_gain = dcgain(minreal(zpk(s*sys)));

% initial conditions
if ss_gain(1,1) == 0; x0_0 = 1; else; x0_0 = ss_ref/ss_gain(1,1); end
x0_0 = x0_0 * x0_ref;
xi_0 = [0 0];

% subsitute A with controlled versions
A0 = A;
A = A_contr;

% impulse response (modal analysis)
if silent < 1
    figure
    impulse(ss(A,B,C,D)), grid on
end



%% Local Controllers
% solve ARE: A'P + PA + Q - PB inv(R) B'P = 0
R = ones(m,m);
Q = eye(n);
P = are(A, B*inv(R)*B', Q);

% compute gain K
K = inv(R)*B'*P;

% compute coupling gain c
in_deg_diag = diag(in_deg);         % in-degrees matrix
L = in_deg_diag - Adj;              % Laplacian matrix of graph
G = diag(g);            % pinning matrix
lambda = eig(L+G);       
c_min = 1 / (2*min(real(lambda)));
c = 2*c_min;



%% Global Control Stability
eig_Ac = [];
for i = 1:N
    eig_Ac = [eig_Ac, eig(A - c*lambda(i)*B*K)'];
end
if silent < 2
    eig_Ac
end



%% Leader Observer
% standard Luenberger observer
L_obs = place(A',C', [-1 -2])';
x0_0_obs = zeros(n,1);

if silent < 2
    eig_Lobs = eig(A-L_obs*C)'
end


%% Local Observers
% observer dynamics
x_hat = sym('x_hat', [n,1], 'real');
ct = sym('ct', [n,1], 'real');          % corrective term
ui = sym('ui', 'real');

x_hat_d = A*x_hat + B*ui - ct;

matlabFunction(x_hat_d, 'File', 'coopObserver', 'Vars', {[x_hat;ui;ct]});
% used in interpred matlab function in simulink 

% inital condition
xi_0_obs = zeros(n,1);

% observer gain
% solve ARE: AP + PA' + Q - PC' inv(R) CP = 0
Ro = ones(m,m);
Qo = eye(n);
Po = are(A', C'*inv(Ro)*C, Qo);

% compute gain F
F = Po*C'*inv(Ro);

% compute coupling gain co
co = c;

% use or not local observer
if local_obs
    F = L_obs / c;      % s.t. c*F = L_obs, standard Lueneberger
end



%% Global Obsever Stability
eig_Ao = [];
for i = 1:N
    if ~local_obs
        eig_Ao = [eig_Ao, eig(A - co*lambda(i)*F*C)'];
    else
        eig_Ao = [eig_Ao, eig(A + co*F*C)'];
    end
end
if silent < 2
    eig_Ao
end


%% Simulation
t_sim = 10;
out = sim('coop_reg_SVFB.slx', 'SrcWorkspace', 'current');



%% Elaborate results

% reshape node states
xi_sim = {N};
for foll_n = 1:N
    xi_sim{foll_n} = out.xi_hat_all(:,(foll_n-1)*n+1:foll_n*n)';
end
% reshape node output errors
y_tilde_hist = {N};
for foll_n = 1:N
    y_tilde_hist{foll_n} = out.yi_tilde_all(:,(foll_n-1)*p+1:foll_n*p)';
end
% collect node outputs
yi_sim = {out.y1(1,:)' out.y2(1,:)' out.y3(1,:)' out.y4(1,:)' out.y5(1,:)' out.y6(1,:)'};

% collect node command inputs
ui_sim = {out.u1(1,:)' out.u2(1,:)' out.u3(1,:)' out.u4(1,:)' out.u5(1,:)' out.u6(1,:)'};

% output values
t_sim = out.tout;
x0_sim = out.x0_hat;
y0_sim = out.y0;

% global disagreement error
delta = zeros(n*N, length(out.tout));
delta_MS = zeros(length(out.tout),1);
for t = 1:length(out.tout)
    x0_bar = kron(ones(N,1), out.x0_hat(t,:)');
    delta(:,t) = out.xi_hat_all(t,:)' - x0_bar;
    delta_MS(t) = 1/N * norm(delta(:,t))^2;
end



%% Plot results
if silent > 0
    return
end

% plot leader dynamics
figure, plot(t_sim, x0_sim), grid on
title('Leader S_0')
legend('x_0','v_0')

% plot nodes states and outputs
for foll_n = 1:N

    figure

    subplot(3,1,1)
    plot(t_sim, xi_sim{foll_n}(1,:)), hold on
    plot(t_sim, x0_sim(:,1), '--'), grid on
    str = sprintf('x_%i', foll_n); legend(str,'x_0')
    str = sprintf('S_%i - position', foll_n); title(str)

    subplot(3,1,2)
    plot(t_sim, xi_sim{foll_n}(2,:)), hold on
    plot(t_sim, x0_sim(:,2), '--'), grid on
    str = sprintf('v_%i', foll_n); legend(str,'v_0')
    str = sprintf('S_%i - velocity', foll_n); title(str)

    subplot(3,1,3)
    plot(t_sim, yi_sim{foll_n}), hold on
    plot(t_sim, y0_sim, '--'), grid on
    str = sprintf('y_%i', foll_n); legend(str,'y_0')
    str = sprintf('S_%i - output', foll_n); title(str)

end

% global disagreement error
figure
plot(out.tout, delta_MS), grid on
title('Mean Square of global disagreement error')