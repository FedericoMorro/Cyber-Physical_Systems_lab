function [y, C, x_hat, eta] = e01_rand_noisy_mes_gen(q, p, k)

% C
C = randn(q,p);     % ~N(0,1)

% x_hat k-sparse
x_hat = zeros(p,1);
while nnz(x_hat) < k        % card(x_hat) < k
    index = randi(p);       % uniformly distributed integer random index
    if x_hat(index) == 0
        uniform_rand = 2*rand(1);   % ~U([0,2])
        if uniform_rand < 1
            x_hat(index) = uniform_rand - 2;    % in [-2,-1]
        else
            x_hat(index) = uniform_rand;        % in [1,2]
        end
    end
end

% noise
sigma = 1e-2;
eta = sigma * randn(q,1);     % ~N(0,sigma^2)

% measurements
y = C*x_hat + eta;

end