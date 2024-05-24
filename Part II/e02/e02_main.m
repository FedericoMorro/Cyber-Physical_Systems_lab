clear
close all
clc

format compact

load exp_data.mat

% parameters
p = 5;
N = size(u,1);
Delta_eta = 5;

% iteration variables
PUI = zeros(5,2);

% iterate over p parameters to estimate PUI_i
for i = 1:p

    % perform computation for min and max of PUI
    for j = 1:2
        
        objPoly.typeCone = 1;
        objPoly.dimVar = N+p;
        objPoly.degree = 1;
        objPoly.noTerms = 1;
        objPoly.supports = zeros(1,N+p);
        objPoly.supports(i) = 1;
        if j == 1, objPoly.coef = 1; else, objPoly.coef = -1; end
        
        for k = 3:N
            ineqPolySys{k-2}.typeCone = -1;
            ineqPolySys{k-2}.dimVar = N+p;
            ineqPolySys{k-2}.degree = 2;
            ineqPolySys{k-2}.noTerms = 9;
            
            supports = zeros(9,N+p);
            supports(2,p+k) = 1;
            supports(3,1) = 1;
            supports(4,1) = 1;  supports(4,p+k-1) = 1;
            supports(5,2) = 1;
            supports(6,2) = 1;  supports(6,p+k-2) = 1;
            supports(7,3) = 1;
            supports(8,4) = 1;
            supports(9,5) = 1;

            ineqPolySys{k-2}.supports = supports;
            ineqPolySys{k-2}.coef = [
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
        
        ubd = [1e10*ones(1,p)   Delta_eta*ones(1,N)];
        lbd = -ubd;
        
        param.relaxOrder = 2;
        param.POPsolver = 'active-set';
        
        [param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = ...
            sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
        PUI(i,j) = POP.objValueL;
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