%% display results
% UNCOMMENT THE BLOCK CORRESPONDING TO THE SIMULATION PERFORMED

f = figure();
f.Position([3 4]) = [1150, 400]; % figure window size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IN FUNCTION OF lambda (tau constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Support recovery rate plot
subplot(1,2,1)
plot(lambda_array, sup_rec_cnt_array, '.--',  'MarkerSize', 20 )
title('Support recovery rate')
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('%')
xlim(lambda_array([1, end]))
ylim([0 100])
ax = gca;
ax.YGrid = 'on';
hold on
plot(lambda_array, mean(sup_rec_cnt_array)*ones(1,N_ROUTINE), '--')

legend('', 'average')

% Convergence time plot (Max, Mean, Min)
subplot(1,2,2)
semilogy(lambda_array, max_conv_time_array, 'd-')
title('Convergence time')
xlabel('$\lambda$', 'Interpreter','latex')
ylabel( '# of iterations' )
xlim(lambda_array([1, end]))
ylim([50 1e5])
ax = gca;
ax.YGrid = 'on';
hold on
semilogy(lambda_array, mean_conv_time_array, 'x-', 'MarkerSize', 10)
semilogy(lambda_array, min_conv_time_array, 'o-')

legend('Max', 'Mean', 'Min')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IN FUNCTION OF tau (tau_lambda constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Support recovery rate plot
subplot(1,2,1)
plot(tau_array, sup_rec_cnt_array, '.--',  'MarkerSize', 20 )
title('Support recovery rate')
xlabel('$\tau$', 'Interpreter','latex')
ylabel('%')
xlim(tau_array([1, end]))
ylim([0 100])
ax = gca;
ax.YGrid = 'on';
hold on
plot(tau_array, mean(sup_rec_cnt_array)*ones(1,N_ROUTINE), '--')

legend('', 'average')

% Convergence time plot (Max, Mean, Min)
subplot(1,2,2)
semilogy(tau_array, max_conv_time_array, 'd-')
title('Convergence time')
xlabel('$\tau$', 'Interpreter','latex')
ylabel( '# of iterations' )
xlim(tau_array([1, end]))
ylim([50 1e5])
ax = gca;
ax.YGrid = 'on';
hold on
semilogy(tau_array, mean_conv_time_array, 'x-', 'MarkerSize', 10)
semilogy(tau_array, min_conv_time_array, 'o-')

legend('Max', 'Mean', 'Min')
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IN FUNCTION OF ROUTINE #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Support recovery rate plot
subplot(1,2,1)
plot(sup_rec_cnt_array, '.--',  'MarkerSize', 20 )
title('Support recovery rate')
xlabel('routine #')
ylabel('%')
xlim([1 10])
xlim(q_array([1, end]))
ylim([0 100])
ax = gca;
ax.YGrid = 'on';
hold on
plot(mean(sup_rec_cnt_array)*ones(1,N_ROUTINE), '--')

legend('', 'average')

% Convergence time plot (Max, Mean, Min)
subplot(1,2,2)
semilogy(max_conv_time_array, 'd-')
title('Convergence time')
xlabel('routine #')

ylabel( '# of iterations' )
xlim([1 10])
xlim(q_array([1, end]))
ylim([50 1e5])
ax = gca;
ax.YGrid = 'on';
hold on
semilogy(mean_conv_time_array, 'x-', 'MarkerSize', 10)
semilogy(min_conv_time_array, 'o-')

legend('Max', 'Mean', 'Min')
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IN FUNCTION OF INCREASING q VALUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Support recovery rate plot
subplot(1,2,1)
plot(q_array, sup_rec_cnt_array, '.--',  'MarkerSize', 20 )
title('Support recovery rate')
xlabel('q')
ylabel('%')
xlim(q_array([1, end]))
ylim([40 100])
ax = gca;
ax.YGrid = 'on';

% Convergence time plot (Max, Mean, Min)
subplot(1,2,2)
semilogy(q_array, max_conv_time_array, 'd-')
title('Convergence time')
xlabel('q')
ylabel( '# of iterations' )
xlim(q_array([1, end]))
ylim([50 1e4])
ax = gca;
ax.YGrid = 'on';
hold on
semilogy(q_array, mean_conv_time_array, 'x-', 'MarkerSize', 10)
semilogy(q_array, min_conv_time_array, 'o-')

legend('Max', 'Mean', 'Min')
%}