function [status, obj, D] = optimize_D_FDMIMO(N, K, H_herm, gamma_min, Pmax, sigma2)

cvx_begin quiet
cvx_solver mosek
variable D(N,K) complex

received_signals  = H_herm * D; % [K, K]
% received_signals  = H_herm * G_diag * A_blk * D;
interference_mat = received_signals .* (1 - eye(K)); % Remove self-interference
% interference_norms = sqrt(sum(abs(interference_mat).^2, 2) + sigma^2);
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

minimize(sum(sum_square_abs(D))*1e6);

subject to
real(diag(received_signals)) >= sqrt(gamma_min) * interference_norms.';
imag(diag(received_signals)) == 0;

cvx_end

if ~strcmp(cvx_status, 'Solved')
    status = false;
    obj = inf;
    D = [];
else
    status = true;
    obj = sum(sum_square_abs(D))*1e6;
end


end

