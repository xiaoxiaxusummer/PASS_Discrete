function [status, obj, sol, ME_gap] = optimize_convex_relax_LMI(M, N, L, K, H_herm, G, G_diag, s_lb, s_ub, gamma_min, Pmax, sigma2)
A_lb_vec = s_lb(1:M); A_lb = reshape(A_lb_vec,[L,N]);
A_ub_vec = s_ub(1:M); A_ub = reshape(A_ub_vec,[L,N]);
U_lb = s_lb(M+1: M+N*K); U_lb=reshape(U_lb,[N,K]);
U_ub = s_ub(M+1: M+N*K);  U_ub = reshape(U_ub,[N,K]);
V_lb = s_lb(M+N*K+1: M+2*N*K); V_lb=reshape(V_lb,[N,K]);
V_ub = s_ub(M+N*K+1: M+2*N*K); V_ub = reshape(V_ub,[N,K]);

cvx_begin quiet
cvx_solver mosek
variable A(L,N) nonnegative
variable D(N,K) complex
variable U(L,L,N)
variable V(K,K,N)
variable mu(M,K) complex

Z = H_herm * G_diag * mu;
interference_mat = Z .* (1 - eye(K)); % Remove self-interference
% interference_norms = sqrt(sum(abs(interference_mat).^2, 2)/sigma^2 + 1);
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

minimize(sum(sum_square_abs(mu))*1e6);
subject to
% sum(sum_square_abs(mu))/L <= Pmax; % (16c)
imag(diag(Z)) == 0; % (17c)
real(diag(Z)) >= sqrt(gamma_min) * interference_norms.'; % (17b)

for n = 1:N
    [U(:,:,n), mu((n-1)*L+1:n*L,:)*1e4, A(:,n); mu((n-1)*L+1:n*L,:)'*1e4, V(:,:,n), D(n,:)'*1e4; A(:,n).', D(n,:)*1e4, 1] == hermitian_semidefinite(L+K+1);
    trace(U(:,:,n) - diag(A(:,n))) <= 0;
    U(:,:,n) == hermitian_semidefinite(L);
    V(:,:,n) == hermitian_semidefinite(K);
end

real(D) >= U_lb; real(D) <= U_ub;
imag(D) >= V_lb; imag(D) <= V_ub;
A <= A_ub; A >= A_lb; % (25e)

cvx_end
ME_gap = +inf;
if ~strcmp(cvx_status, 'Solved')
    status = false;
    obj = inf;
    sol = struct();
else
    status = true;
    % D = U + 1j* V;
    obj = sum(sum(abs(mu).^2))/L;
    sol = struct('Z',Z, 'D',D, 'A', A, 'mu', mu, 'A_diag', diag(reshape(A,[M,1])));
end


end