function [A_des_eig, x0_ref, t_fin] = reference_config(ref_type, magn, freq)
% desired eigenvalues, reference and simulation duration

A_des_eig = 0.5*[0 -1];         % constant
x0_ref = [magn 0]';
t_fin = 10;

if ref_type == "ramp"
    A_des_eig = [0 0];          % ramp
    x0_ref = [0 magn]';
    t_fin = 15;
elseif ref_type == "sinusoidal"
    A_des_eig = [0+freq*1i 0-freq*1i];  % sinusoidal
    t_fin = 20;
end

end