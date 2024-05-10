% run e05_main.m to set up workspace variables

close all
clc

format compact

f = figure(1);
f.Position([3 4]) = [1150, 400];

subplot(1,2,1)
grid on, hold on
title('avg{ 1/2 ||x_{est} - x||_1 }'), xlabel('k'), ylabel('err_x')

subplot(1,2,2)
grid on, hold on
title('avg{ 1/2 ||a_{est} - a||_1 }'), xlabel('k'), ylabel('err_a')

k_init = 2;
max_k_x = round(max(k_x_conv) * 1.1);
max_k_a = round(max(k_a_conv) * 1.1);

% plot evolution of error in norm
for j = 1:length(Q_vec)
    
    subplot(1,2,1)
    plot(k_init:max_k_x, 1/2*x_norm_error(j,k_init:max_k_x))
    
    subplot(1,2,2)
    plot(k_init:max_k_a, 1/2*a_norm_error(j,k_init:max_k_a))

end

subplot(1,2,1)
legend(Q_names)
subplot(1,2,2)
legend(Q_names)