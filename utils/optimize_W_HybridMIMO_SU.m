function W = optimize_W_HybridMIMO_SU(K, M, F, D, H_herm, rho, lambda_W, gamma_min, sigma2)
cvx_begin quiet
cvx_solver mosek
variable W(M,K) complex

received_signals  = H_herm * W; % a complex number

residual_W = W - F*D;
minimize(sum(sum_square_abs(W)) + 1/2/rho*sum(sum_square_abs(residual_W+rho*lambda_W)));
subject to
real(received_signals) >= sqrt(gamma_min) * sqrt(sigma2);
imag(received_signals) == 0;
cvx_end


end