clear
close all
clc

format compact

% load results
load sim_res.mat

names = {'RMS{\epsilon}','t_{ga}','||u||_2^2','t_{s,2}','\delta_y'};

figure
for i = 1:5
    subplot(1,5,i), hold on
    for ref_n = 1:length(refs)
        plot(c_fact_vec, c_fact_d{ref_n}(i,:))
    end
    xlim([max(c_fact_vec)*0.95, max(c_fact_vec)*1.05]) , grid on, title(names{i})
end
legend('step','ramp','sin')