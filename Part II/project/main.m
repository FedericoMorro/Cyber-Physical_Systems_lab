clear all
close all
clc

format compact

% silent: 0 plot and print, 1 print and plot external, 2 print external
%   -1 for (also) impulse response
silent = 2;

% reference type
refs = {"step", "ramp", "sinusoidal"};

% local or global observers
local_obs = false;

% algorithm parameters
n = 2; m = 1;
par.c_fact = 2;
par.co_fact = 1;
par.R = eye(m);
par.Q = eye(n);
par.Ro = eye(m);
par.Qo = eye(n);

% parameters vectors
N_SIM = 10;
c_fact_vec = linspace(1,10,N_SIM);
QR_vec = logspace(-2,2,N_SIM);
co_fact_vec = logspace(0,6,N_SIM);
QoRo_vec = logspace(-6,6,N_SIM);



%% c_fact, distributed observers

c_fact_d = {length(refs)};
for ref_n = 1:length(refs)
    c_fact_d{ref_n} = [];

    for curr_c = c_fact_vec
        par.c_fact = curr_c;
        ref_type = refs{ref_n};

        metrics = CPS_sim(ref_type, par, false, silent);

        c_fact_d{ref_n} = [c_fact_d{ref_n} metrics];
    end
end

par.c_fact = 2;



%% Q R, distributed observers

QR_d = {length(refs)};
for ref_n = 1:length(refs)
    QR_d{ref_n} = [];

    for curr_QR = QR_vec
        par.Q = curr_QR * eye(n);
        ref_type = refs{ref_n};

        metrics = CPS_sim(ref_type, par, false, silent);

        QR_d{ref_n} = [QR_d{ref_n} metrics];
    end
end

par.Q = eye(n);



%% co_fact, distributed observers

co_fact = {length(refs)};
for ref_n = 1:length(refs)
    co_fact{ref_n} = [];

    for curr_co = co_fact_vec
        par.co_fact = curr_co;
        ref_type = refs{ref_n};

        metrics = CPS_sim(ref_type, par, false, silent);

        co_fact{ref_n} = [co_fact{ref_n} metrics];
    end
end

par.co_fact = 1;



%% Qo Ro, distributed observers

QoRo = {length(refs)};
for ref_n = 1:length(refs)
    QoRo{ref_n} = [];

    for curr_QoRo = QoRo_vec
        par.Qo = curr_QoRo * eye(n);
        ref_type = refs{ref_n};

        metrics = CPS_sim(ref_type, par, false, silent);

        QoRo{ref_n} = [QoRo{ref_n} metrics];
    end
end

par.Qo = eye(n);



%% c_fact, local observers

c_fact_l = {length(refs)};
for ref_n = 1:length(refs)
    c_fact_l{ref_n} = [];

    for curr_c = c_fact_vec
        par.c_fact = curr_c;
        ref_type = refs{ref_n};

        metrics = CPS_sim(ref_type, par, true, silent);

        c_fact_l{ref_n} = [c_fact_l{ref_n} metrics];
    end
end

par.c_fact = 2;



%% Q R, local observers

QR_l = {length(refs)};
for ref_n = 1:length(refs)
    QR_l{ref_n} = [];

    for curr_QR = QR_vec
        par.Q = curr_QR * eye(n);
        ref_type = refs{ref_n};

        metrics = CPS_sim(ref_type, par, true, silent);

        QR_l{ref_n} = [QR_l{ref_n} metrics];
    end
end

par.Q = eye(n);



%% Save
save sim_res