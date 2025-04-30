function [obj, sol] = optimize_d_SU(L, N, K, A, Q, h_herm, G, G_diag, gamma_min, sigma2, L_sel)
M = N*L;
P_opt = gamma_min * sigma2 / sum(abs(h_herm * diag(A(:)) * G).^2 ./ L_sel); 
d_opt = sqrt(P_opt) * (h_herm * diag(A(:)) * G)' / sqrt(sum(abs(h_herm * diag(A(:)) * G).^2));
obj = watts_to_dbm(P_opt);
% obj = P_opt * 1e8;
sol = struct('A', A, 'd', d_opt, 'P', P_opt);
end
