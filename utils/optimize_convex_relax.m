function [status, obj, sol] = optimize_convex_relax(M, N, L, K, H_herm, G, G_diag, s_lb, s_ub, gamma_min, sigma2)
A_lb_vec = s_lb(1:M); A_lb = reshape(A_lb_vec,[L,N]);
A_ub_vec = s_ub(1:M); A_ub = reshape(A_ub_vec,[L,N]);
U_lb = s_lb(M+1: M+N*K); U_lb=reshape(U_lb,[N,K]);
U_ub = s_ub(M+1: M+N*K);  U_ub = reshape(U_ub,[N,K]);
V_lb = s_lb(M+N*K+1: M+2*N*K); V_lb=reshape(V_lb,[N,K]);
V_ub = s_ub(M+N*K+1: M+2*N*K); V_ub = reshape(V_ub,[N,K]);

cvx_begin quiet
cvx_solver mosek
variable A(L,N) nonnegative
variable U(N,K)
variable V(N,K)
variable mu(M,K) complex


received_signals  = H_herm * G_diag * mu; % [K, K]
interference_mat = received_signals .* (1 - eye(K)); % Remove self-interference
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

minimize(sum(sum_square_abs(mu))*1e7);
subject to
real(diag(received_signals)) >= sqrt(gamma_min) * interference_norms.';
imag(diag(received_signals)) == 0;

for n = 1: N
    real(mu((n-1)*L+1:n*L,:))* 1e6 >= (A(:,n)*U_lb(n,:) + A_lb(:,n)*U(n,:) - A_lb(:,n)*U_lb(n,:))* 1e6;
    real(mu((n-1)*L+1:n*L,:))* 1e6 >= (A(:,n)*U_ub(n,:) + A_ub(:,n)*U(n,:) - A_ub(:,n)*U_ub(n,:))* 1e6;
    real(mu((n-1)*L+1:n*L,:))* 1e6 <= (A(:,n)*U_ub(n,:) + A_lb(:,n)*U(n,:) - A_lb(:,n)*U_ub(n,:))* 1e6;
    real(mu((n-1)*L+1:n*L,:))* 1e6 <= (A(:,n)*U_lb(n,:) + A_ub(:,n)*U(n,:) - A_ub(:,n)*U_lb(n,:))* 1e6;
    imag(mu((n-1)*L+1:n*L,:))* 1e6 >= (A(:,n)*V_lb(n,:) + A_lb(:,n)*V(n,:) - A_lb(:,n)*V_lb(n,:))* 1e6;
    imag(mu((n-1)*L+1:n*L,:))* 1e6 >= (A(:,n)*V_ub(n,:) + A_ub(:,n)*V(n,:) - A_ub(:,n)*V_ub(n,:))* 1e6;
    imag(mu((n-1)*L+1:n*L,:))* 1e6 <= (A(:,n)*V_ub(n,:) + A_lb(:,n)*V(n,:) - A_lb(:,n)*V_ub(n,:))* 1e6;
    imag(mu((n-1)*L+1:n*L,:))* 1e6 <= (A(:,n)*V_lb(n,:) + A_ub(:,n)*V(n,:) - A_ub(:,n)*V_lb(n,:))* 1e6;
end

U >= U_lb; U <= U_ub;
V >= V_lb; V <= V_ub;
A <= A_ub; A >= A_lb;

cvx_end

if ~strcmp(cvx_status, 'Solved')
    status = false;
    obj = inf;
    sol = struct();
else
    status = true;
    D = U + 1j* V;
    obj = watts_to_dbm(sum(sum(abs(mu).^2)));
    sol = struct('D',D, 'A', A, 'mu', mu, 'A_diag', diag(reshape(A,[M,1])));
end


end