function [A_des_eig] = reference_config(ref_type, freq)
% desired eigenvalues

A_des_eig = 0.5*[0 -1];         % constant

if ref_type == "ramp"
    A_des_eig = [0 0];          % ramp
elseif ref_type == "sinusoidal"
    A_des_eig = [0+freq*1i 0-freq*1i];  % sinusoidal
end

end