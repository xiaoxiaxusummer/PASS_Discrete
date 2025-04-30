function [status, obj, sol] = optimize_D_argument(L, N, K, H_herm, G, G_diag, A, gamma_min, Pmax, sigma2, s_lb, s_ub)
M = N*L;
A_lb = s_lb(1:M); 
A_ub = s_ub(1:M);
P_lb = s_lb(M+1: M+N*K); P_lb = reshape(P_lb,[N,K]);
P_ub = s_ub(M+1: M+N*K); P_ub = reshape(P_ub,[N,K]);
theta_lb = reshape(s_lb(M+N*K+1: M+N*K*2),[N,K]);
theta_ub = reshape(s_ub(M+N*K+1: M+N*K*2),[N,K]);

U = (cos(theta_lb) + cos(theta_ub))/2; 
V = (sin(theta_lb) + sin(theta_ub))/2;
S = sqrt(U.^2 + V.^2);

for n = 1:N
    idx = find((theta_ub(n,:) - theta_lb(n,:)) > pi);
    U(n,idx) = 0; V(n,idx) = 0; S(n,idx) = 0;
end

cvx_begin quiet
cvx_solver mosek
variable D(N,K) complex
variable P(N,K) nonnegative

expression Z(M,K)
for n = 1: N
    Z((n-1)*L+1:n*L,1:K) = A(:,n)*D(n,:);
end
received_signals  = H_herm * G_diag * Z; % [K, K]
interference_mat = received_signals .* (1 - eye(K)); % Remove self-interference
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

minimize(sum(sum_square_abs(Z))*1e6);

subject to
real(diag(received_signals)) >= sqrt(gamma_min) * interference_norms.';
imag(diag(received_signals)) == 0;
(U .* real(D) + V .* imag(D)) * 1e3 >= P .* S * 1e3;
real(D).^2 + imag(D).^2 <= P;
P >= P_lb; P <= P_ub;

cvx_end

P = real(D).^2 + imag(D).^2;

if ~strcmp(cvx_status, 'Solved')
    status = false;
    obj = inf;
    sol = struct();
else
    status = true;
    obj = sum(sum_square_abs(full(Z)))*1e6;
    sol = struct('P', P, 'A', A, 'Z', Z, 'D', D, 'A_diag', diag(reshape(A,[M,1])));
end


end

