function [D_opt, obj] = MIMO_SU(K, N, H_herm, gamma_min, sigma2)
P = sigma2 * gamma_min / sum(abs(H_herm).^2);
D_opt = sqrt(P) * H_herm' / sqrt(sum(abs(H_herm).^2));  
obj = watts_to_dbm(P);
end