clear
close all
clc

format compact

% paramters
p = 100;    % # of cells
q = 25;     % # of sensors
load("localization.mat")        % A, D, y


% estimated vectors
x = z(1:n);
a = z(n+1:n+q);