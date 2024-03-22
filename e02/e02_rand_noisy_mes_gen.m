function [y, C, x_hat, a, eta] = e02_rand_noisy_mes_gen(n, q, h, aware)

% C
C = randn(q,n);     % ~N(0,1)

% x_hat
x_hat = randn(n,1);     % ~N(0,1)

% uncorrupted measurements
y_ideal = C*x_hat;

% noise
sigma = 1e-2;
eta = sigma * randn(q,1);     % ~N(0,sigma^2)

% a
a = zeros(q,1);
while nnz(a) < h        % card(a) < h
    index = randi(q);       % uniformly distributed integer random index
    if a(index) == 0
        % different generation for un/aware attacks
        if aware
            a(index) = y_ideal(index) / 2;
        else
            uniform_rand = 2*rand(1);   % ~U([0,2])
            if uniform_rand < 1
                a(index) = uniform_rand - 2;    % in [-2,-1]
            else
                a(index) = uniform_rand;        % in [1,2]
            end
        end

    end
end

% measurements
y = y_ideal + eta + a;

end