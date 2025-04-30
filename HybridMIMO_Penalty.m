function [F, D, obj] =  HybridMIMO_Penalty(K, M, N, L, H_herm, gamma_min, sigma2, subconnected)

rho = 0.05;
F = randn(M,N) * 1j + randn(M,N);
F = F./abs(F);
D = randn(N,K) * 1j + randn(N,K);
W = randn(M,K) * 1j + randn(M,K);
lambda_W = zeros(M,K) * 1j + zeros(M,K);

Gamma = kron(eye(N), ones(L,1));
F = F .* Gamma; 
f_all = [];  residuals = [];

Tmax = 30;
BCD_itermax = 20; 
BCD_tol = 1e-4;
for t = 1:Tmax
    for i = 1:BCD_itermax
        Q = H_herm * F * D;
        received_signals = sum(abs(Q).^2,1) + sigma2;
        effective_response = conj(diag(Q));
        numerators = abs(effective_response).^2;
        SINRs = numerators./(received_signals-numerators);

        W = optimize_W_HybridMIMO(K, M, F, D, H_herm, lambda_W, rho, gamma_min, sigma2);
        F = optimize_F_HybridMIMO(K, N, M, W, F, D, lambda_W, rho, subconnected);

        D = optimize_D_HybridMIMO(K, N, W, F, lambda_W, rho, H_herm, sigma2, gamma_min);
        

        f_all(1+length(f_all)) = watts_to_dbm(sum(sum_square_abs(F*D)));
        residual_W = W - F*D;
        
        if i >= 3 && abs(f_all(end)-f_all(end-1)) <= BCD_tol
            break;
        end
    end
    % disp(f_all(end));
    if t > 1
        last_max_residual = max_residual;
    end
    max_residual = max(sum(abs(residual_W),1));
    residuals(1+length(residuals)) = max_residual;
    if t>=3 && (abs(f_all(end) - f_all(end-1)) < 1e-3 || max_residual<1e-5)
        break;
    end
end

obj = f_all(end);

end