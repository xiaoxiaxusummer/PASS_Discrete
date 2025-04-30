function [power_dbm,sol,power_log,residuals] = run_continuous_search(M, N, L, K, loc_U, A, phi, cfg, discrete_A, discrete_x, discrete_W)

alg_param = struct();
alg_param.rho = 0.1;
alg_param.Tmax = 300;

gamma_min = cfg.gamma_min;
sigma2 = cfg.sigma2;
n_eff = cfg.n_eff;
beta = cfg.beta;
kappa = cfg.kappa;
Delta = cfg.lambda/2;
h_PAA = cfg.h_PAA;
S_PAA = cfg.S_x;
Gamma = zeros(M,N);
for n = 1:N
    Gamma((n-1)*L+1:n*L,n) = discrete_A(:,n);
end
Gamma = Gamma ./ sqrt(sum(discrete_A,1));
x = reshape(discrete_x, [M,1]);
x_init = x;
R = sqrt((x - A).^2 + phi.^2); 
Theta = kappa * (R+ n_eff .* x);
U = sqrt(beta) ./ R .* exp(-1j*Theta);

W = discrete_W;
Q =  U.' * Gamma * W;
interference_mat = Q .* (1 - eye(K)); % Remove self-interference
SINRs = abs(diag(Q)).^2 ./ (sum(abs(interference_mat).^2,2) + sigma2);
power_dbm = watts_to_dbm(sum(sum_square_abs(W)));
disp(['initial: ', num2str(power_dbm)]);

rho = alg_param.rho; 
power_log = [power_dbm]; residuals = [];
for t = 1:alg_param.Tmax
    %  Start searching
    num_particles = 50;
    for n = 1:N
        for l = 1:L
            x_mat = reshape(x,[L,N]);
            if l > 1
                x_lb = max(x_mat(l-1,n) - Delta,0);
            else
                x_lb = 0;
            end
            if l < L
                x_ub = min(x_mat(l+1,n) - Delta,S_PAA);
            else
                x_ub = S_PAA;
            end
            loc_particles = linspace(x_lb,x_ub,num_particles);
            for p = 1:num_particles
                x_new = x_mat; x_new(l,n) = loc_particles(p); x_new = reshape(x_new,[M,1]);
                R = sqrt((x_new - A).^2 + phi.^2); 
                Theta = kappa * (R+ n_eff .* x_new);
                U = sqrt(beta) ./ R .* exp(-1j*Theta);
                [W_new, status] = optimize_D_continuous(K, N, U, Gamma, sigma2, gamma_min);
                power_new = watts_to_dbm(sum(sum_square_abs(W_new)));
                if status && power_new < power_dbm
                    x = x_new; x_mat = reshape(x_new,[L,N]); W = W_new;
                    power_dbm = power_new;
                    % disp(['update: ', num2str(power_dbm)]);
                end
            end
        end
    end

    acc_power_iter = watts_to_dbm(sum(sum_square_abs(W)));
    power_log = [power_log, acc_power_iter];
    % disp(['acc power:',num2str(acc_power_iter)]);
    if t>=2 && abs(power_log(end) - power_log(end-1)) <= 1e-3
        break;
    end
end
disp(['gamma_min: ', num2str(cfg.gamma_min_dB), ', K: ', num2str(K), ', L: ',num2str(L), ...
    ', space: ', num2str(cfg.S_x), ', final power: ', num2str(acc_power_iter)]);

sol = struct( 'x', x,  'W', W,  'U', U,  'Theta', Theta, ...
        'loc_user', loc_U, 'accurate_power', acc_power_iter,  'rate_log', power_log);
end