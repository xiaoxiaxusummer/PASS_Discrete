function [status, obj, sol] = optimize_relax_SU(L, N, K, Q, h_herm, G, gamma_min, sigma2, s_lb, s_ub, L_sel)
M = N*L;
a_lb = s_lb(1:M); A_lb = reshape(a_lb,[L,N]);
a_ub = s_ub(1:M); A_ub = reshape(a_ub,[L,N]);


cvx_begin quiet
cvx_solver mosek
cvx_precision best
variable Z(L,L,N) nonnegative symmetric
variable A(L,N) nonnegative
variable eta(1,N) nonnegative

maximize(sum(eta./L_sel)/gamma_min/sigma2*1e7);
subject to
    for n = 1:N
        Z(:,:,n) >= A(:,n) * A_lb(:,n).' + A_lb(:,n) * A(:,n).' - A_lb(:,n) * A_lb(:,n).';
        Z(:,:,n) >= A(:,n) * A_ub(:,n).' + A_ub(:,n) * A(:,n).' - A_ub(:,n) * A_ub(:,n).';
        Z(:,:,n) <= A(:,n) * A_lb(:,n).' + A_ub(:,n) * A(:,n).' - A_ub(:,n) * A_lb(:,n).';
        Z(:,:,n) <= A(:,n) * A_ub(:,n).' + A_lb(:,n) * A(:,n).' - A_lb(:,n) * A_ub(:,n).';
        real(Q(:,n)' * Z(:,:,n) * Q(:,n)) == eta(n);
    end
    A >= A_lb; A <= A_ub;
    sum(A,1) == L_sel;

cvx_end


if strcmp(cvx_status, 'Solved') 
    status = true;
    P = gamma_min * sigma2 / sum(eta ./ L_sel);
    d_opt = sqrt(P) * (h_herm * diag(A(:)) * G)' / sqrt(sum(abs(h_herm * diag(A(:)) * G).^2));
    obj = watts_to_dbm(P);
    % obj = P * 1e8;
    sol = struct('P', P, 'A', A, 'Z', Z, 'd', d_opt);
elseif strcmp(cvx_status, 'Inaccurate/Solved')
    status = true;
    P = gamma_min * sigma2  / sum(eta ./ L_sel);
    d_opt = sqrt(P) * (h_herm * diag(A(:)) * G)' / sqrt(sum(abs(h_herm * diag(A(:)) * G).^2));
    obj = watts_to_dbm(P);
    % obj = P * 1e8;
    sol = struct('P', P, 'A', A, 'Z', Z, 'd', d_opt);
else
    status = false;
    obj = inf;
    sol = struct();
end


end

