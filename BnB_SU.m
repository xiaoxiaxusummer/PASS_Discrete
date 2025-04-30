function [A_opt, d_opt, GUB, GLB_hist, GUB_hist] = BnB_SU(K, M, N, L, h_herm, G, G_diag, gamma_min, sigma2, L_sel)
G_mat = reshape(diag(G_diag),[L,N]);
H_mat = conj(reshape(h_herm(:),[L,N]));
Q = zeros(L,N) * 1j;
for n = 1:N
    Q(:,n) = diag(G_mat(:,n)') * H_mat(:,n);
end

% Initial relaxation: Solve MINLP as NLP by relaxing x \in {0,1} to [0,1]
A_lb = zeros(M,1);  % Lower bound for binary variables
A_ub = ones(M,1);    % Upper bound for binary variables

s_lb_all = [A_lb]; % [M+2NK+2K, num_box]
s_ub_all = [A_ub]; % [M+2NK+2K, num_box]

% Initial solution
[status, f_lb, sol_relax] = optimize_relax_SU(L, N, K, Q, h_herm, G, gamma_min, sigma2, s_lb_all, s_ub_all, L_sel);
if status
    GLB = f_lb;
    GLB_hist = [GLB];
    f_lb_all = [f_lb]; % local upper bound
    sol_relax_all = {sol_relax};
else
    disp('Infeasible problem');
end

% Find a feasible solution
[sorted_values, sorted_indices] = sort(sol_relax.A, 1, 'descend'); % Sort each column in descending order
A_proj = zeros(size(sol_relax.A)); % Create a zero matrix of same size
for n = 1:N
    A_proj(sorted_indices(1:L_sel(n), n), n) = 1; % Set the top L indices to 1
end
[f_ub, sol_proj] = optimize_d_SU(L, N, K, A_proj, Q, h_herm, G, G_diag, gamma_min, sigma2, L_sel);
GUB = f_ub;
GUB_hist = [GUB];
f_ub_all = [f_ub]; % local upper bound;
sol_all = {sol_proj};


while ~isempty(s_lb_all)
    % select the best node from candidate list
    [~, idx] = min(f_lb_all);
    s_lb = s_lb_all(:,idx); s_ub = s_ub_all(:,idx); sol_relax = sol_relax_all{idx}; sol_proj = sol_all{idx};
    s_lb_all(:,idx) = []; s_ub_all(:,idx) = [];
    f_lb_all(idx) = []; f_ub_all(idx) = []; sol_relax_all(idx) = []; sol_all(idx) = [];

    % branching
    s_lb_left = s_lb; s_ub_left = s_ub; s_lb_right = s_lb; s_ub_right = s_ub;
    % [~, e] = max([s_ub(1:M) - s_lb(1:M)]);
    if full(sum(abs(sol_relax.A(:) - round(sol_relax.A(:)))))>0.1
        [~, e] = min(abs(sol_relax.A(:)-0.5));
    else
        [~, e] = max([s_ub(1:M) - s_lb(1:M)]);
    end
    s_lb_left(e) = 0; s_ub_left(e) = 0;  s_lb_right(e) = 1; s_ub_right(e) = 1;
    
    s_lb_child = [s_lb_left, s_lb_right]; s_ub_child = [s_ub_left, s_ub_right];

    % Solve relaxed problem
    for child = 1:2
        if any(sum(reshape(s_lb_child(:,child),[L,N]),1) > L_sel) || any(sum(reshape(s_ub_child(:,child),[L,N]),1) < L_sel)
            continue;
        end
        [status, f_lb, sol_relax] = optimize_relax_SU(L, N, K, Q, h_herm, G_diag, gamma_min, sigma2, s_lb_child(:,child), s_ub_child(:,child), L_sel);
        % If infeasible or non-optimal, prune
        if ~status || f_lb - 1e-2 > GUB 
            continue;
        end

        % Find a feasible solution
        [sorted_values, sorted_indices] = sort(sol_relax.A, 1, 'descend'); % Sort each column in descending order
        A_proj = zeros(size(sol_relax.A)); % Create a zero matrix of same size
        for n = 1:N
            A_proj(sorted_indices(1:L_sel(n), n), n) = 1; % Set the top L indices to 1
        end
        [f_ub, sol_proj] = optimize_d_SU(L, N, K, A_proj, Q, h_herm, G, G_diag, gamma_min, sigma2, L_sel);
        if f_lb - 1e-2 > f_ub
            disp("upper bound warning");
        end
        % Update the currently best solution
        if f_ub <= GUB
            GUB = f_ub;
            A_opt = sol_proj.A; d_opt = sol_proj.d;
            pruned_idx = any(f_lb_all > GUB, 1);
            f_lb_all(pruned_idx) = []; f_ub_all(pruned_idx) = []; sol_relax_all(pruned_idx) = []; sol_all(pruned_idx) = [];
            s_lb_all(:,pruned_idx) = []; s_ub_all(:,pruned_idx) = [];
        end

        % prune if fathomed
        if abs(f_lb - f_ub) <= 1e-4
            continue;
        end
        s_lb_all = [s_lb_all, s_lb_child(:,child)];
        s_ub_all = [s_ub_all, s_ub_child(:,child)];
        f_lb_all = [f_lb_all, f_lb]; f_ub_all = [f_ub_all, f_ub];
        sol_relax_all{length(sol_relax_all)+1} = sol_relax; sol_all{length(sol_all)+1} = sol_proj;
    end
    GLB = min(f_lb_all); 
    GLB_hist = [GLB_hist, GLB]; 
    GUB_hist = [GUB_hist, GUB];
    % disp(['GUB:', num2str(GUB), '; GLB:', num2str(GLB), '; number of remaining boxes:', num2str(length(f_lb_all))]);
end
end
