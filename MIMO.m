function [D_opt, obj] = MIMO(K, N, H_herm, Pmax, gamma_min, sigma2)
[status, obj, D_opt] = optimize_D_FDMIMO(N, K, H_herm, gamma_min, Pmax, sigma2);
obj = watts_to_dbm(obj/1e6);
end