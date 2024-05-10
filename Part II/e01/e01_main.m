clear
close all
clc

format compact


s = tf('s');

% plant specification
Gp = 100/(s^2 + 1.2*s + 1);

% discretized model
dt = 1;                     % sample time
Gd = c2d(Gp, dt, 'zoh');
Dd = tf(1,Gd.Denominator{1},dt);
% Gd(z) = (th_3 z^2 + th_4 z + th_5) / (z^2 + th_1 z + th_2)
%   z -> q  = (th_3 q^2 + th_4 q + th_5) / (q^2 + th_1 q + th_2)
%           = (th_3 + th_4 q^-1 + th_5 q^-2) / (1 + th_1 q^-1 + th_2 q^-2)
%           = N(q^-1) / D(q^-1)

% extract exact values of parameters
theta_true = [Gd.Denominator{1}(2:end) Gd.Numerator{1}]';

% 2nd order system -> mininimum H = 5, since 5 parameters to be estimated
H_min = 7;
H_max = 100;

% error variance
sigma = 5;



%% Noiseless LS
% iterations variables
theta_nn_err = zeros(H_max-H_min+1, 1);

for H = H_min:H_max

    % random input generation
    u = rand(H,1);
    
    % simulated output 
    y = lsim(Gd,u);

    % solve LS
    A = [-y(2:H-1) -y(1:H-2) u(3:H) u(2:H-1) u(1:H-2)];
    b = y(3:H);
    theta_nn = pinv(A)*b;

    % calculate error in norm
    theta_nn_err(H-H_min+1) = norm(theta_nn - theta_true, 2)^2;

end

% theta_nn_err(H) = 0,  forall H    -> no plot



%% LS with equation error e
% follow LS assumption

% y(k) = Gd(q^-1)*u(k)
% y_tilde(k)*D(q^-1) = N(q^-1)*u(k) + e(k)
% y_tilde(k) = N(q^-1)/D(q^-1) * u(k) + 1/D(q^-1) * e(k)
% strange: LS assumption states that uncertainty affecting the output is
%   adding a term which depends on the denominator of the system which is
%   unknown since we are trying to estimate it

% iterations variables
theta_ee_err = zeros(H_max-H_min+1, 1);

for H = H_min:H_max

    % random input generation
    u = rand(H,1);

    % simulated output with gaussian equation error
    e = sigma * randn(H,1);
    y = lsim(Gd,u);
    y_tilde = y + lsim(Dd,e);

    % solve LS
    A = [-y_tilde(2:H-1) -y_tilde(1:H-2) u(3:H) u(2:H-1) u(1:H-2)];
    b = y_tilde(3:H);
    theta_ee = pinv(A)*b;

    % calculate error in norm
    theta_ee_err(H-H_min+1) = norm(theta_ee - theta_true, 2)^2;

end

% plot
figure
plot(H_min:H_max, theta_ee_err)
grid on
hold on


%% LS with output measurement error
% correct real experiment setup

% y(k) = Gd(q^-1)*u(k)
% y_tilde(k) = y(k) + eta(k) = N(q^-1)/D(q^-1) * u(k) + eta(k)

% iterations variables
theta_oe_err = zeros(H_max-H_min+1, 1);

for H = H_min:H_max

    % random input generation
    u = rand(H,1);

    % simulated output with gaussian equation error
    eta = sigma * randn(H,1);
    y = lsim(Gd,u);
    y_tilde = y + eta;

    % solve LS
    A = [-y_tilde(2:H-1) -y_tilde(1:H-2) u(3:H) u(2:H-1) u(1:H-2)];
    b = y_tilde(3:H);
    theta_oe = pinv(A)*b;

    % calculate error in norm
    theta_oe_err(H-H_min+1) = norm(theta_oe - theta_true, 2)^2;

end

% plot
plot(H_min:H_max, theta_oe_err)
grid on