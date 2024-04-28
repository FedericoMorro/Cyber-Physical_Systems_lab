%% display results

f = figure();
f.Position([3 4]) = [1200, 400]; % figure window size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IN FUNCTION OF ROUTINE #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Support recovery rate plot
subplot(1,2,1)
plot(att_detection_array, '.--',  'MarkerSize', 20 )
title('Attack detection rate')
xlabel('routine #')
ylabel('%')
xlim([1 N_ROUTINE])
ylim([0 100])
ax = gca;
ax.YGrid = 'on';
hold on
plot(mean(att_detection_array)*ones(1,N_ROUTINE), '--')

legend('', 'average')

% Convergence time plot (Max, Mean, Min)
subplot(1,2,2)
plot(max_estim_acc_array, 'd-')
title('Estimation accuracy')
xlabel('routine #')

ylabel( '$\| \tilde{x} - \hat{x} \|_2^2$', Interpreter='latex' )
xlim([1 N_ROUTINE])
%ylim([50 1e5])
ax = gca;
ax.YGrid = 'on';
hold on
plot(mean_estim_acc_array, 'x-', 'MarkerSize', 10)
plot(min_estim_acc_array, 'o-')

legend('Max', 'Mean', 'Min')
