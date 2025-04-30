function [status, obj, sol] = optimize_convex_relax2(M, N, L, K, H_herm, G, G_diag, s_lb, s_ub, gamma_min, Pmax, sigma2)
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

U = repelem(U, L, 1);
V = repelem(V, L ,1);
S = repelem(S, L, 1);

cvx_begin quiet
cvx_solver mosek
variable A(L,N) nonnegative
variable P(N,K) nonnegative
variable Z(M,K) complex
variable R(M,K) nonnegative

received_signals  = H_herm * G_diag * Z; % [K, K]
interference_mat = received_signals .* (1 - eye(K)); % Remove self-interference
% interference_norms = sqrt(sum(abs(interference_mat).^2, 2)/sigma^2 + 1);
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

minimize(sum(sum(R))*1e6);
% minimize(sum(sum_square_abs(Z))*1e6);
subject to
real(diag(received_signals)) >= sqrt(gamma_min) * interference_norms.'; % (17b)
imag(diag(received_signals)) == 0; % (17c)
(real(Z).^2 + imag(Z).^2) * 1e5 <= R * 1e5;

(U .* real(Z) + V .* imag(Z)) * 1e3 >= R .* S * 1e3;


for n = 1: N
    R((n-1)*L+1:n*L,:) * 1e5 >= (A(:,n) * P_lb(n,:)) * 1e5 ;
    R((n-1)*L+1:n*L,:) * 1e5 >= (A(:,n) * P_ub(n,:) + ones(L,1) * P(n,:) - ones(L,1) * P_ub(n,:)) * 1e5;
    R((n-1)*L+1:n*L,:) * 1e5 <= (A(:,n) * P_ub(n,:)) * 1e5;
    R((n-1)*L+1:n*L,:) * 1e5 <= (A(:,n) * P_lb(n,:) + ones(L,1) * P(n,:) - ones(L,1) * P_lb(n,:)) * 1e5;
end

P >= P_lb; P <= P_ub;
A(:) <= A_ub; A(:) >= A_lb; 
% sum(A,1) == round(L*0.5);

cvx_end

if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
    status = false;
    obj = inf;
    sol = struct();
elseif strcmp(cvx_status, 'Solved')
    status = true;
    obj = sum(sum(R))*1e6;
    sol = struct('P', P, 'A', A, 'Z', Z, 'R', R, 'A_diag', diag(reshape(A,[M,1])));
else
    status = true;
    obj = sum(sum(R))*1e6;
    sol = struct('P', P, 'A', A, 'Z', Z, 'R', R, 'A_diag', diag(reshape(A,[M,1])));
end


end