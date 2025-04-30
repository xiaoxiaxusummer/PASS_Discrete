function [f_obj, residual_U, residual_Q, residual_Theta] = print_log(t, i, Q, v_E, alpha, U, D, Theta, x, R, Gamma, rho,  ... 
                                                       lambda_U, lambda_Q, lambda_Theta, rate, sigma2, n_eff, beta, kappa, A, phi, h_PAA)

    MSE = cal_MSE_by_Q(Q, v_E, sigma2);
    WMMSE = sum(alpha.*MSE-log2(alpha));
    residual_U = U.*R - sqrt(beta) .*exp(-1j*Theta);
    residual_Q = Q - U.' * Gamma * D;
    residual_Theta = Theta-  kappa*R - kappa*n_eff * x;
    L_U = sum(sum(abs(residual_U + rho*lambda_U).^2));
    L_Q = sum(sum(abs(residual_Q + rho*lambda_Q).^2));
    L_Theta = sum(sum(abs(residual_Theta + rho*lambda_Theta).^2));
    f_obj = WMMSE + 1/(2*rho)*(L_U+L_Q+L_Theta);
       
    fprintf('Out iter: %d, BCD iter: %d; WMMSE: %.12f; f_obj: %.12f; R_U: %.12f; R_Q: %.12f; R_Theta: %.12f; sum rate: %.12f\n', ...
                t, i, WMMSE, f_obj, sum(abs(residual_U(:))), sum(abs(residual_Q(:))), sum(abs(residual_Theta(:))), sum(rate(:)));

    % 调用函数 cal_accurate_rate 并打印结果
    acc_Rate = cal_accurate_rate(x, D, beta, n_eff, kappa, sigma2, Gamma, A, phi, h_PAA);
    fprintf('accurate sum rate: %.12f\n', sum(acc_Rate(:)));

end

function MSE = cal_MSE_by_Q(Q, v_E, sigma2)
    signal = abs(v_E).^2 .* sum(abs(Q).^2,2);
    noise_term = sigma2 * (abs(v_E).^2);
    cross_term = -2 * real(v_E .* diag(Q));
    MSE = signal + noise_term + 1 + cross_term;
end