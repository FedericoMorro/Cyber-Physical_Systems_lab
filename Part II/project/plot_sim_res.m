clear
close all
clc

format compact

% load results
load sim_res.mat

names = {'RMS({\epsilon})','t_{ga}','||u||_2^2','t_{s, 2%}','\delta_y','x_{bar}'};
y = {c_fact_d, QR_d, co_fact, QoRo, c_fact_l, QR_l};
x = {c_fact_vec, QR_vec, co_fact_vec, QoRo_vec, c_fact_vec, QR_vec};
x_names = {'c_f', 'Q/R', 'co_f', 'Qo/Ro', 'c_f', 'Q/R'};


for ind = 1:6
    fig_n = ind;
    if ind == 5; fig_n = 1; elseif ind == 6; fig_n = 2; end
    style = '-';
    if ind == 5 || ind == 6; style = '--'; end

    f = figure(fig_n);
    f.Position([3,4]) = [1200, 600];
    for i = 1:6
        subplot(2,3,i)
        for ref_n = 1:length(refs)
            if ind == 1 || ind == 5
                plot(x{ind}, y{ind}{ref_n}(i,:),style)
            else
                semilogx(x{ind}, y{ind}{ref_n}(i,:),style)
            end
            hold on
        end
        axis('padded'), grid on
        title(names{i})
        if i > 3
            xlabel(x_names{ind})
        end
    end

    if ~(ind == 5 || ind == 6)
        legend('step','ramp','sin')
    else
        legend('step','ramp','sin', 'step loc','ramp loc','sin loc')
    end
end