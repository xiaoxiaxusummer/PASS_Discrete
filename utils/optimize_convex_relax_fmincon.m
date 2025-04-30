function [status, obj, sol] = optimize_convex_relax_fmincon(M, N, L, K, H_herm, G, G_diag, s_lb, s_ub, gamma_min, Pmax, sigma2)

    % **1. 变量界限解析**
    A_lb_vec = s_lb(1:M); A_lb = reshape(A_lb_vec,[L,N]);
    A_ub_vec = s_ub(1:M); A_ub = reshape(A_ub_vec,[L,N]);
    U_lb = reshape(s_lb(M+1:M+N*K), [N,K]);
    U_ub = reshape(s_ub(M+1:M+N*K), [N,K]);
    V_lb = reshape(s_lb(M+N*K+1:M+2*N*K), [N,K]);
    V_ub = reshape(s_ub(M+N*K+1:M+2*N*K), [N,K]);

    % **2. 变量向量化**
    num_A = L * N;
    num_U = N * K;
    num_V = N * K;
    num_mu = M * K;
    num_vars = num_A + num_U + num_V+ num_mu*2;  % 总优化变量数

    % **3. 初始点**
    A0 = (A_ub+A_lb)/2;
    % D = randn(N,K) + randn(N,K)*1j;
    equiv_channel = (H_herm * diag(A0(:)) * G)'; % [N,K]
    lambda_coeff = diag(ones(K,1)*Pmax/K); % [K,K]
    D = (eye(N) + 1/sigma2 * equiv_channel * lambda_coeff * equiv_channel')\equiv_channel;
    power_all = sqrt(ones(1,K)*Pmax/K) ./ sum(abs(D),1); % (1,K)
    D = D .* power_all;
    U = real(D); V = imag(D); 
    mu_real = zeros(M,K); mu_imag = zeros(M,K);
    for n = 1:N
        idx = (n-1)*L+1 : n*L;
        mu_real(idx,:) = A0(:,n) * U(n,:);
        mu_imag(idx,:) = A0(:,n) * V(n,:);
    end
    x0 = [(A_ub_vec+A_lb_vec)/2; U(:); V(:); mu_real(:); mu_imag(:)];

    % **4. 变量边界**
    lb = [A_lb(:); U_lb(:); V_lb(:); -inf(num_mu*2,1)];
    ub = [A_ub(:); U_ub(:); V_ub(:);  inf(num_mu*2,1)];
    
    % **5. 定义目标函数**
    objective = @(x) obj_func(x, M, N, L, K);

    % **6. 非线性约束**
    nonlcon  = @(x) nonlcon_func(x, M, N, L, K, H_herm, G_diag, gamma_min, Pmax, sigma2, A_lb, A_ub, U_lb, U_ub, V_lb, V_ub);

    % **8️. 调用 fmincon**
% options = optimoptions('fmincon', ...
%     'Display', 'iter', ...
%     'Algorithm', 'interior-point', ...
%     'MaxFunctionEvaluations', 1e5, ...
%     'MaxIterations', 1e4, ...
%     'OptimalityTolerance', 1e-6, ...
%     'StepTolerance', 1e-10, ...
%     'ConstraintTolerance', 1e-6, ...
%     'HonorBounds', true, ...             % 确保变量边界严格满足
%     'ScaleProblem', true, ...            % **自动缩放问题，提高数值稳定性**
%     'SubproblemAlgorithm', 'factorization', ... % **改进子问题求解**
%     'EnableFeasibilityMode', true);      % **启用可行性模式**

    options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...              % ✅ 使用 SQP 代替 interior-point
    'MaxFunctionEvaluations', 1e5, ...
    'MaxIterations', 1e4, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-6);


    [x_opt, fval, exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, options);

    % **9. 结果处理**
    if exitflag <= 0
        status = false;
        obj = inf;
        sol = struct();
        return;
    end

    % 解析变量
    A = reshape(x(1:num_A), [L, N]);
    U = reshape(x(num_A+1:num_A+num_U), [N, K]);
    V = reshape(x(num_A+num_U+1:num_A+num_U+num_V), [N, K]);
    mu_real = reshape(x(num_A+num_U+num_V+1:num_A+num_U+num_V+num_mu), [M, K]);
    mu_imag = reshape(x(num_A+num_U+num_V+num_mu+1:end), [M, K]);
    mu = mu_real + 1j*mu_imag;

    D = U + 1j * V;
    Z = H_herm * G_diag * mu;
    obj = sum(sum(abs(D).^2));

    sol = struct('Z', Z, 'D', D, 'A', A, 'mu', mu, 'A_diag', diag(reshape(A, [M,1])));
    status = true;

end


function f = obj_func(x, M, N, L, K)
    num_A = L * N;
    num_U = N * K;
    num_V = N * K;

    % 提取 U 和 V
    U = reshape(x(num_A+1:num_A+num_U), [N, K]);
    V = reshape(x(num_A+num_U+1:num_A+num_U+num_V), [N, K]);

    % 目标函数：最小化功率
    f = sum(sum(abs(U).^2)) + sum(sum(abs(V).^2));
end


function [c, ceq] = nonlcon_func(x, M, N, L, K, H_herm, G_diag, gamma_min, Pmax, sigma2, A_lb, A_ub, U_lb, U_ub, V_lb, V_ub)

    num_A = L * N;
    num_mu = M * K;
    num_U = N * K;
    num_V = N*K;

    % 解析变量
    A = reshape(x(1:num_A), [L, N]);
    U = reshape(x(num_A+1:num_A+num_U), [N, K]);
    V = reshape(x(num_A+num_U+1:num_A+num_U+num_V), [N, K]);
    mu_real = reshape(x(num_A+num_U+num_V+1:num_A+num_U+num_V+num_mu), [M, K]);
    mu_imag = reshape(x(num_A+num_U+num_V+num_mu+1:end), [M, K]);


    % 复数重构
    mu = mu_real + 1j * mu_imag;
    Z = H_herm * G_diag * mu;

    % **强制 `Z` 对角线为实数**
    ceq = imag(diag(Z));  % (17c) 保证 Z 的对角线是实数

    % 计算干扰项
    interference_mat = Z .* (1 - eye(K));
    interference_norms = sqrt(sum(abs(interference_mat).^2, 2) + sigma2);

    % 约束条件
    c = [
        sum(sum(abs(U).^2)) + sum(sum(abs(V).^2)) - Pmax;  % (16c) 总功率约束
        sqrt(gamma_min) * interference_norms - real(diag(Z)); % (17b) SINR 约束
    ];

    for n = 1:N
        idx = (n-1)*L+1 : n*L;

        % 处理实部 McCormick 约束
        c = [c;
            reshape(mu_real(idx,:) - (A(:,n) * U_lb(n,:) + A_lb(:,n) * U(n,:) - A_lb(:,n) * U_lb(n,:)), [], 1);
            reshape(mu_real(idx,:) - (A(:,n) * U_ub(n,:) + A_ub(:,n) * U(n,:) - A_ub(:,n) * U_ub(n,:)),[],1);
            reshape(-(mu_real(idx,:) - (A(:,n) * U_ub(n,:) + A_lb(:,n) * U(n,:) - A_lb(:,n) * U_ub(n,:))),[],1);
            reshape(-(mu_real(idx,:) - (A(:,n) * U_lb(n,:) + A_ub(:,n) * U(n,:) - A_ub(:,n) * U_lb(n,:))),[],1);
        ];

        % 处理虚部 McCormick 约束
        c = [c;
            reshape(mu_imag(idx,:) - (A(:,n) * V_lb(n,:) + A_lb(:,n) * V(n,:) - A_lb(:,n) * V_lb(n,:)),[],1);
            reshape(mu_imag(idx,:) - (A(:,n) * V_ub(n,:) + A_ub(:,n) * V(n,:) - A_ub(:,n) * V_ub(n,:)),[],1);
            reshape(-(mu_imag(idx,:) - (A(:,n) * V_ub(n,:) + A_lb(:,n) * V(n,:) - A_lb(:,n) * V_ub(n,:))),[],1);
            reshape(-(mu_imag(idx,:) - (A(:,n) * V_lb(n,:) + A_ub(:,n) * V(n,:) - A_ub(:,n) * V_lb(n,:))),[],1);
        ];
    end

end
