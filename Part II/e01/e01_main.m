clear
close all
clc

format compact

s = tf('s');

% plant specification
Gp = 100/(s^2 + 1.2*s + 1);

% discretized model
dt = 1;                % sample time
Gd = c2d(Gp, dt, 'zoh');

% random input generation
H = 50;
u = rand(H, 1);

% simulated output 
t_sim = dt * 1:H;
w = lsim(Gd, u, t_sim);
