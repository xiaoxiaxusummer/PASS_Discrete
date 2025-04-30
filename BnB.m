function [A_opt, D_opt, GUB, GLB_hist, GUB_hist] = BnB(K, M, N, L, H_herm, G, G_diag, P0, gamma_min, sigma2)
A0 = ones(L,N);  A_diag = diag(reshape(A0,[M,1]));

% Initial relaxation: Solve MINLP as NLP by relaxing x \in {0,1} to [0,1]
A_lb = zeros(M,1);  % Lower bound for binary variables
A_ub = ones(M,1);    % Upper bound for binary variables
U_lb = ones(N*K,1) * (-sqrt(P0)); V_lb = ones(N*K,1) * (-sqrt(P0));
U_ub = ones(N*K,1) * sqrt(P0); V_ub = ones(N*K,1) * sqrt(P0);

s_lb_all = [A_lb; U_lb; V_lb]; % [M+2NK+2K, num_box]
s_ub_all = [A_ub; U_ub; V_ub]; % [M+2NK+2K, num_box]

% Initial solution
[status, f_ub, D_proj] = optimize_D(L, N, K, H_herm, G_diag, A0, gamma_min, sigma2);
if status
    GUB = f_ub;
    GUB_hist = [GUB];
    A_opt = A0;
    D_opt = D_proj;
    f_ub_all = [f_ub]; % local upper bound
end

[status, f_lb, sol] = optimize_convex_relax(M, N, L, K, H_herm, G, G_diag, s_lb_all, s_ub_all, gamma_min, sigma2);
if status
    GLB = f_lb;
    GLB_hist = [GLB];
    sol_GLB = sol;
    f_lb_all = [f_lb]; % local upper bound
else
    disp('Infeasible problem');
end

while ~isempty(s_lb_all)
    % [~, idx] = min(f_ub_all - f_lb_all); % select the best node from candidate list
    [~, idx] = min(f_lb_all); % select the best node from candidate list
    s_lb = s_lb_all(:,idx); s_ub = s_ub_all(:,idx);
    s_lb_all(:,idx) = []; s_ub_all(:,idx) = [];
    f_lb_all(idx) = []; f_ub_all(idx) = [];

    % edge selection
    [~, e] = max([s_ub(1:M) - s_lb(1:M); (s_ub(M+1:end) - s_lb(M+1:end))./(sqrt(P0)*2)]);
    % branching
    s_lb_left = s_lb; s_ub_left = s_ub; s_lb_right = s_lb; s_ub_right = s_ub;
    if e <= M
        s_lb_left(e) = 0; s_ub_left(e) = 0;  s_lb_right(e) = 1; s_ub_right(e) = 1;
    else
        s_ub_left(e) = (s_lb(e)+s_ub(e))/2;
        s_lb_right(e) = (s_lb(e)+s_ub(e))/2;optimize_D
    end
    s_lb_child = [s_lb_left, s_lb_right]; s_ub_child = [s_ub_left, s_ub_right];

    % Solve relaxed problem
    for child = 1:2
        [status, f_lb, sol] = optimize_convex_relax(M, N, L, K, H_herm, G, G_diag, s_lb_child(:,child), s_ub_child(:,child), gamma_min, sigma2);
        % If infeasible or non-optimal, prune
        if ~status || f_lb - 1e-4  > GUB
            continue;
        end

        % Find a feasible solution
        % A_proj = round(sol.A);
        A_proj = reshape(s_ub_child(1:M,child), [L, N]);
        [status, f_ub, D_proj] = optimize_D(L, N, K, H_herm, G_diag, A_proj, gamma_min, sigma2);

        % Update the currently best solution
        if f_ub <= GUB
            GUB = f_ub;
            A_opt = A_proj;
            D_opt = D_proj;
            pruned_idx = any(f_lb_all > GUB, 1);
            f_lb_all(pruned_idx) = []; f_ub_all(pruned_idx) = [];
            s_lb_all(:,pruned_idx) = []; s_ub_all(:,pruned_idx) = [];
        end

        % prune if fathomed
        if abs(f_lb - f_ub) <= 1e-4
            continue;
        end
        s_lb_all = [s_lb_all, s_lb_child(:,child)];
        s_ub_all = [s_ub_all, s_ub_child(:,child)];
        f_lb_all = [f_lb_all, f_lb]; f_ub_all = [f_ub_all, f_ub];
    end
    GLB = min(f_lb_all);
    GLB_hist = [GLB_hist, GLB]; GUB_hist = [GUB_hist, GUB];
    % disp(['GUB:', num2str(GUB), '; GLB:', num2str(GLB), '; number of remaining boxes:', num2str(length(f_lb_all))]);

end
end
