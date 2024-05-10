clear all
clc
close all

format compact

n = 1000;
k_fin = 70;

res_mse = zeros(k_fin, 3);
res_conv = zeros(2, 3);
res_a_mse = zeros(k_fin, 3);

thr = 1.2;

tic;
x_mse_avg = gpuArray(zeros(k_fin, 1));
a_mse_avg = gpuArray(zeros(k_fin, 1));
x_exact_time = 0;
x_avg_conv = 0;
parfor i = 1:n
    [x_err, fc, a_err] = rand_atk(2, k_fin, thr);
    x_mse_avg = x_mse_avg + x_err;
    a_mse_avg = a_mse_avg + a_err;

    if fc
        x_exact_time = x_exact_time + 1;
        x_avg_conv = x_avg_conv + fc;
    end
end
x_mse_avg = x_mse_avg/n;
x_avg_conv = x_avg_conv/x_exact_time;
a_mse_avg = a_mse_avg/n;
res_mse(:,1) = gather(x_mse_avg);
res_conv(1, 1) = x_avg_conv;
res_conv(2, 1) = x_exact_time;
res_a_mse(:, 1) = gather(a_mse_avg);
t1 = toc

tic;
x_mse_avg = gpuArray(zeros(k_fin, 1));
a_mse_avg = gpuArray(zeros(k_fin, 1));
x_exact_time = 0;
x_avg_conv = 0;
parfor i = 1:n
    [x_err, fc, a_err] = rand_atk(4, k_fin, thr);
    x_mse_avg = x_mse_avg + x_err;
    a_mse_avg = a_mse_avg + a_err;
    if fc
        x_exact_time = x_exact_time + 1;
        x_avg_conv = x_avg_conv + fc;
    end
end
x_mse_avg = x_mse_avg/n;
x_avg_conv = x_avg_conv/x_exact_time;
a_mse_avg = a_mse_avg/n;
res_mse(:, 2) = gather(x_mse_avg);
res_conv(1, 2) = x_avg_conv;
res_conv(2, 2) = x_exact_time;
res_a_mse(:, 2) = gather(a_mse_avg);
t2 = toc

tic;
x_mse_avg = gpuArray(zeros(k_fin, 1));
a_mse_avg = gpuArray(zeros(k_fin, 1));
x_exact_time = 0;
x_avg_conv = 0;
parfor i = 1:n
    [x_err, fc, a_err] = rand_atk(6, k_fin, thr);
    x_mse_avg = x_mse_avg + x_err;
    a_mse_avg = a_mse_avg + a_err;
    if fc
        x_exact_time = x_exact_time + 1;
        x_avg_conv = x_avg_conv + fc;
    end
end
x_mse_avg = x_mse_avg/n;
x_avg_conv = x_avg_conv/x_exact_time;
a_mse_avg = a_mse_avg/n;
res_mse(:, 3) = gather(x_mse_avg);
res_conv(1, 3) = x_avg_conv;
res_conv(2, 3) = x_exact_time;
res_a_mse(:, 3) = gather(a_mse_avg);
t3 = toc

save("test_report_const_atks");