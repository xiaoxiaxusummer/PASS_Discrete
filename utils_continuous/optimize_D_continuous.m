function [D, status] = optimize_D_continuous(K, N, U, Gamma, sigma2, gamma_min)
status = true;
cvx_begin quiet
cvx_solver mosek
    variable D(N, K) complex
    
    received_signals  =  U.' * Gamma * D; % [K, K]

    interference_mat = received_signals .* (1 - eye(K)); % Remove self-interference
    interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
    interference_norms = norms(interference_norms.',2);  % [1,K]
    
    minimize(sum(sum_square_abs(D))*1e6);
    subject to
        real(diag(received_signals)) >= sqrt(gamma_min) * interference_norms.';
        imag(diag(received_signals)) == 0;
        
cvx_end

if ~strcmp(cvx_status, 'Solved')
    status = false;
end
end