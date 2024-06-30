clear
close all
clc

format compact

% load results
load sim_res.mat

names = {'RMS{\epsilon}','t_{ga}','||u||_2^2','t_{s,2}','\delta_y','x_{bar}'};
y = {c_fact_d, QR_d, co_fact, QoRo, c_fact_l, QR_l};
x = {c_fact_vec, QR_vec, co_fact_vec, QoRo_vec, c_fact_vec, QR_vec};
x_names = {'c_f', 'Q/R', 'co_f', 'Qo/Ro', 'c_f', 'Q/R'};


for ind = 1:6
    f = figure();
    f.Position([3,4]) = [1200, 600];
    for i = 1:6
        subplot(2,3,i)
        for ref_n = 1:length(refs)
            if ind == 1 || ind == 5
                plot(x{ind}, y{ind}{ref_n}(i,:))
            else
                semilogx(x{ind}, y{ind}{ref_n}(i,:))
            end
            hold on
        end
        axis('padded'), grid on
        title(names{i}), xlabel(x_names{ind})
    end
    legend('step','ramp','sin')
end