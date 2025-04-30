function [status, obj, sol] = optimize_relax_SU_SW(L, N, K, h_herm, G, G_diag, gamma_min, Pmax, sigma2, s_lb, s_ub, L_sel)
M = N*L;
a_lb = s_lb(1:L); a_lb = reshape(a_lb,[L,1]);
a_ub = s_ub(1:L); a_ub = reshape(a_ub,[L,1]);

equiv_channel = h_herm * G_diag;


cvx_begin quiet
cvx_solver sdpt3
cvx_precision best
% variable Z(M,M) symmetric
variable Q(L,L) nonnegative symmetric
variable a(L,1) nonnegative
variable eta nonnegative

maximize(eta/gamma_min/sigma2/L_sel*1e6);
% minimize(gamma_min*sigma2*sum(pow_p(eta,-1))*1e8);
subject to
    real(equiv_channel * Q * equiv_channel') == eta;
    Q >= a * a_lb.' + a_lb * a.' - a_lb * a_lb.'; 
    Q >= a * a_ub.' + a_ub * a.' - a_ub * a_ub.';
    Q <= a * a_lb.' + a_ub * a.' - a_lb * a_ub.';
    Q <= a * a_ub.' + a_lb * a.' - a_ub * a_lb.'; 
    a >= a_lb; a <= a_ub;
    sum(a) == L_sel;

cvx_end


if strcmp(cvx_status, 'Solved') 
    status = true;
    P = gamma_min * sigma2 * L_sel / eta;
    d_opt = sqrt(P) * (h_herm * G_diag * a)' / abs(h_herm * G_diag * a);
    % obj = watts_to_dbm(P);
    obj = P * 1e8;
    sol = struct('P', P, 'a', a, 'Z', Z, 'd', d_opt);
elseif strcmp(cvx_status, 'Inaccurate/Solved')
    status = true;
    P = gamma_min * sigma2 * L_sel / eta;
    d_opt = sqrt(P) * (h_herm * G_diag * a)' / abs(h_herm * G_diag * a);
    % obj = watts_to_dbm(P);
    obj = P * 1e8;
    sol = struct('P', P, 'a', a, 'Z', Z, 'd', d_opt);
else
    status = false;
    obj = inf;
    sol = struct();
end


end

