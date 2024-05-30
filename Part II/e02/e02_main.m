clear
close all
clc

format compact

load exp_data.mat

% parameters
n = 2;      % system order
p = 5;      % number of parameters
N = size(u,1);
Delta_eta = 5;
par_bound = 30;

% iteration variables
PUI = zeros(5,2);

% iterate over p parameters to estimate PUI_i
for i = 1:p

    % perform computation for min and max of PUI
    for j = 1:2
        
        objPoly.typeCone = 1;
        objPoly.dimVar = N+p;   % var = [th(1) ... th(p) eta(1) ... eta(N)]
        objPoly.degree = 1;
        objPoly.noTerms = 1;
        objPoly.supports = zeros(1,N+p);
        objPoly.supports(i) = 1;
        if j == 1, objPoly.coef = 1; else, objPoly.coef = -1; end
        
        for k = n+1:N
            ineqPolySys{k-n}.typeCone = -1;
            ineqPolySys{k-n}.dimVar = N+p;
            ineqPolySys{k-n}.degree = 2;
            ineqPolySys{k-n}.noTerms = 9;
            
            supports = zeros(9,N+p);
            supports(2,p+k) = 1;
            supports(3,1) = 1;
            supports(4,1) = 1;  supports(4,p+k-1) = 1;
            supports(5,2) = 1;
            supports(6,2) = 1;  supports(6,p+k-2) = 1;
            supports(7,3) = 1;
            supports(8,4) = 1;
            supports(9,5) = 1;

            ineqPolySys{k-n}.supports = supports;
            ineqPolySys{k-n}.coef = [
                y_tilde(k)
                -1
                y_tilde(k-1)
                -1
                y_tilde(k-2)
                -1
                -u(k)
                -u(k-1)
                -u(k-2)
            ];
        end
        
        ubd = [par_bound*ones(1,p)   Delta_eta*ones(1,N)];
        lbd = -ubd;
        
        param.relaxOrder = 2;
        param.POPsolver = 'active-set';
        
        [param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = ...
            sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);

        if j == 1
            PUI(i,j) = POP.objValueL;
        else
            PUI(i,j) = - POP.objValueL;
        end
        %POP.xVectL

    end

end
PUI

theta_c = zeros(p,1);
for i = 1:p
    theta_c(i) = mean(PUI(i,:));
end
theta_c

% to force dc-gain equal to 118, add another constraint
dc_gain = 118;