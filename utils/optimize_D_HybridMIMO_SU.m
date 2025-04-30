function d = optimize_D_HybridMIMO_SU(K, N, W, F, lambda_W, rho, H_herm, sigma2, gamma_min)
cvx_begin quiet
cvx_solver mosek
variable d(N,1) complex

residual_W = W - F*d;

received_signals  = H_herm * F * d; % a complex number

minimize(1/2/rho*sum(sum_square_abs(residual_W+lambda_W*rho)));


subject to
real(received_signals) >= sqrt(gamma_min) * sqrt(sigma2);
imag(received_signals) == 0;

cvx_end

end