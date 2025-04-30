function WP_list = preference_list(L, N, K, H_herm, G, G_diag, sigma2, gamma_min, A_new)
    [status, P_new, D, W] = optimize_D_matching(L, N, K, H_herm, G, G_diag, A_new, gamma_min, Pmax, sigma2);


end

