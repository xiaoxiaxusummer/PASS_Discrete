function [utility, D, W] = optimize_D_matching(L, N, K, H_herm, G, G_diag, A, gamma_min, sigma2)
M = N*L;

cvx_begin quiet
cvx_solver mosek
variable D(N,K) complex

expression mu(M,K)
for n = 1: N
    mu((n-1)*L+1:n*L,1:K) = A(:,n)*D(n,:);
end
received_signals  = H_herm * G_diag * mu; % [K, K]
interference_mat = received_signals .* (1 - eye(K)); % Remove self-interference
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

minimize(sum(sum_square_abs(mu))*1e7);

subject to
real(diag(received_signals)) >= sqrt(gamma_min) * interference_norms.';
imag(diag(received_signals)) == 0;

cvx_end

if ~strcmp(cvx_status, 'Solved')
    status = false;
    utility = inf;
    D = []; W = [];
else
    status = true;
    % utility = sum(sum_square_abs(full(mu)))*1e6;
    utility = watts_to_dbm(sum(sum_square_abs(full(mu))));
    W = D.*sqrt(sum(A,1).');
end


end

