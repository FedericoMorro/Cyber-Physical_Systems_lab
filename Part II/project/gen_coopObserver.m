function x_hat_d = gen_coopObserver(in1)
%gen_coopObserver
%    X_HAT_D = gen_coopObserver(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    30-Jun-2024 11:40:49

ct1 = in1(4,:);
ct2 = in1(5,:);
ui = in1(3,:);
x_hat1 = in1(1,:);
x_hat2 = in1(2,:);
x_hat_d = [-ct1+x_hat2;-ct2-ui.*9.9453-x_hat1./4.0+x_hat2.*8.763486191350272e-17];
end
