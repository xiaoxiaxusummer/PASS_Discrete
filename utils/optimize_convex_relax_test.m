function [status, obj, sol, ME_gap] = optimize_convex_relax_test(M, N, L, K, H_herm, G, G_diag, s_lb, s_ub, gamma_min, Pmax, sigma2)
A_lb_vec = s_lb(1:M); A_lb = reshape(A_lb_vec,[L,N]);
A_ub_vec = s_ub(1:M); A_ub = reshape(A_ub_vec,[L,N]);
U_lb = s_lb(M+1: M+N*K); U_lb=reshape(U_lb,[N,K]);
U_ub = s_ub(M+1: M+N*K);  U_ub = reshape(U_ub,[N,K]);
V_lb = s_lb(M+N*K+1: M+2*N*K); V_lb=reshape(V_lb,[N,K]);
V_ub = s_ub(M+N*K+1: M+2*N*K); V_ub = reshape(V_ub,[N,K]);

% blkA_lb = zeros(M,N); blkA_ub = zeros(M,N);
% for n = 1:N
%     blkA_lb((n-1)*L+1:n*L,n) = A_lb(:,n);
%     blkA_ub((n-1)*L+1:n*L,n) = A_ub(:,n);
% end

cvx_begin quiet
cvx_solver mosek
variable A(L,N) nonnegative
variable U(N,K)
variable V(N,K)
variable mu(M,K) complex
% expression blkA(M,N)

Z = H_herm * G_diag * mu;
interference_mat = Z .* (1 - eye(K)); % Remove self-interference
% interference_norms = sqrt(sum(abs(interference_mat).^2, 2)/sigma^2 + 1);
interference_norms = cat(2, interference_mat, ones(K,1)*sqrt(sigma2)); % [K,K+1]
interference_norms = norms(interference_norms.',2);  % [1,K]

% for n = 1:N
%     blkA((n-1)*L+1:n*L,n) = A(:,n);
%     if n > 1
%         blkA(1:(n-1)*L,n) = 0;
%     end
%     if n < N
%         blkA(n*L+1:M,n) = 0;
%     end
% end

minimize(sum(sum_square_abs(U)+sum_square_abs(V)));

subject to
sum(sum_square_abs(U)+sum_square_abs(V)) <= Pmax; % (16c)
imag(diag(Z)) == 0; % (17c)
real(diag(Z)) >= sqrt(gamma_min) * interference_norms.'; % (17b)

% real(mu) >= blkA * U_lb + blkA_lb * U + blkA_lb * U_lb;
% real(mu) >= blkA * U_ub + blkA_ub * U + blkA_ub * U_ub;
% real(mu) <= blkA * U_ub + blkA_lb * U + blkA_lb * U_ub;
% real(mu) <= blkA * U_lb + blkA_ub * U + blkA_ub * U_lb;
% 
% 
% imag(mu) >= blkA * V_lb + blkA_lb * V + blkA_lb * V_lb;
% imag(mu) >= blkA * V_ub + blkA_ub * V + blkA_ub * V_ub;
% imag(mu) <= blkA * V_ub + blkA_lb * V + blkA_lb * V_ub;
% imag(mu) <= blkA * V_lb + blkA_ub * V + blkA_ub * V_lb;


% for n = 1: N
%     for k = 1:K
%         real(mu((n-1)*L+1:n*L,k)) >= A(:,n)*U_lb(n,k) + A_lb(:,n)*U(n,k) - A_lb(:,n)*U_lb(n,k);
%         real(mu((n-1)*L+1:n*L,k)) >= A(:,n)*U_ub(n,k) + A_ub(:,n)*U(n,k) - A_ub(:,n)*U_ub(n,k);
%         real(mu((n-1)*L+1:n*L,k)) <= A(:,n)*U_ub(n,k) + A_lb(:,n)*U(n,k) - A_lb(:,n)*U_ub(n,k);
%         real(mu((n-1)*L+1:n*L,k)) >= A(:,n)*U_lb(n,k) + A_ub(:,n)*U(n,k) - A_ub(:,n)*U_lb(n,k);
%
%         imag(mu((n-1)*L+1:n*L,k)) >= A(:,n)*V_lb(n,k) + A_lb(:,n)*V(n,k) - A_lb(:,n)*V_lb(n,k);
%         imag(mu((n-1)*L+1:n*L,k)) >= A(:,n)*V_ub(n,k) + A_ub(:,n)*V(n,k) - A_ub(:,n)*V_ub(n,k);
%         imag(mu((n-1)*L+1:n*L,k)) <= A(:,n)*V_ub(n,k) + A_lb(:,n)*V(n,k) - A_lb(:,n)*V_ub(n,k);
%         imag(mu((n-1)*L+1:n*L,k)) >= A(:,n)*V_lb(n,k) + A_ub(:,n)*V(n,k) - A_ub(:,n)*V_lb(n,k);
%     end
% end

for n = 1: N
    real(mu((n-1)*L+1:n*L,:)) >= A(:,n)*U_lb(n,:) + A_lb(:,n)*U(n,:) - A_lb(:,n)*U_lb(n,:);
    real(mu((n-1)*L+1:n*L,:)) >= A(:,n)*U_ub(n,:) + A_ub(:,n)*U(n,:) - A_ub(:,n)*U_ub(n,:);
    real(mu((n-1)*L+1:n*L,:)) <= A(:,n)*U_ub(n,:) + A_lb(:,n)*U(n,:) - A_lb(:,n)*U_ub(n,:);
    real(mu((n-1)*L+1:n*L,:)) <= A(:,n)*U_lb(n,:) + A_ub(:,n)*U(n,:) - A_ub(:,n)*U_lb(n,:);

    imag(mu((n-1)*L+1:n*L,:)) >= A(:,n)*V_lb(n,:) + A_lb(:,n)*V(n,:) - A_lb(:,n)*V_lb(n,:);
    imag(mu((n-1)*L+1:n*L,:)) >= A(:,n)*V_ub(n,:) + A_ub(:,n)*V(n,:) - A_ub(:,n)*V_ub(n,:);
    imag(mu((n-1)*L+1:n*L,:)) <= A(:,n)*V_ub(n,:) + A_lb(:,n)*V(n,:) - A_lb(:,n)*V_ub(n,:);
    imag(mu((n-1)*L+1:n*L,:)) <= A(:,n)*V_lb(n,:) + A_ub(:,n)*V(n,:) - A_ub(:,n)*V_lb(n,:);

end

U >= U_lb; U <= U_ub;
V >= V_lb; V <= V_ub;
A <= A_ub; A >= A_lb; % (25e)

cvx_end

if ~strcmp(cvx_status, 'Solved')
    status = false;
    obj = inf;
    sol = struct();
else
    status = true;
    D = U + 1j* V;
    obj = sum(sum(abs(D).^2));
    sol = struct('Z',Z, 'D',D, 'A', A, 'mu', mu, 'A_diag', diag(reshape(A,[M,1])));
    % D_lb = U_lb+1j*V_lb; D_ub = U_ub+1j*V_ub;
    % mu_real_min = max(blkA .* U_lb + blkA_lb .* U + blkA_lb .* U_lb, blkA .* U_ub + blkA_ub .* U + blkA_ub .* U_lb);
    % mu_real_max = min(blkA .* U_lb + blkA_ub .* U + blkA_ub .* U_lb, blkA .* U_ub + blkA_lb .* U + blkA_lb .* U_ub);
    % mu_imag_min = max(blkA .* V_lb + blkA_lb .* V + blkA_lb .* V_lb, blkA .* V_ub + blkA_ub .* V + blkA_ub .* V_ub);
    % mu_imag_max = min(blkA .* V_ub + blkA_lb .* V + blkA_lb .* V_ub, blkA .* V_lb + blkA_ub .* V + blkA_ub .* V_lb);
    % ME_gap = [mu_real_max - mu_real_min; mu_imag_max-mu_imag_min]; % [2M, K]
end


end