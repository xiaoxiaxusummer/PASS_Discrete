function W = optimize_W_HybridMIMO(K, M, F, D, H_herm, lambda_W, rho, gamma_min, sigma2)
cvx_begin quiet
cvx_solver mosek
variable W(M,K) complex

received_signals  = H_herm * W; % [K, K]
interference_mat = received_signals .* (1 - eye(K)); % Remove self-interference
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

residual_W = W - F*D;
minimize(sum(sum_square_abs(W)) + 1/2/rho*sum(sum_square_abs(residual_W+lambda_W*rho)));
subject to
real(diag(received_signals)) >= sqrt(gamma_min) * interference_norms.';
imag(diag(received_signals)) == 0;
cvx_end


end