function [A_opt, D_opt, cost] =  Matching(K, M, N, L, H_herm, G, G_diag, gamma_min, sigma2)
gain_mean = reshape(mean(abs(H_herm).^2,1).',[L,N]);
A0 = zeros(L,N); 
for n = 1:N
    [~,ind]=max(gain_mean(:,n));
    A0(ind,n) = 1; 
end
A_opt = A0;
[cost, D_opt, W_opt] = optimize_D_matching(L, N, K, H_herm, G, G_diag, A_opt, gamma_min, sigma2);

cost = [cost]; iter = 0;
% Welfare-driven matching
while true
    cost_t = cost(end); iter = iter + 1;
    for n = 1:N
        for l= 1:L
            if A_opt(l,n) == 0
                % From (l,0)&(n,0) -> (l,n)&(0,0)
                A_new = A_opt; A_new(l,n) = 1;
                [cost_new, D, W] = optimize_D_matching(L, N, K, H_herm, G, G_diag, A_new, gamma_min, sigma2);
                if cost_new < cost
                    A_opt = A_new; D_opt = D; W_opt = W; cost=[cost, cost_new];
                end
            else
                %%
                % case 1: delete a matched pair: (l,n)&(0,0) -> (l',n)&(l,0)
                % case 2: replace a matched pair: (l,n)&(l',0) -> (l',n)&(l,0)
                % case 3: swap matched pairs: (l,n)&(l',n') -> (l,n')&(l',n)
                %%
                for n1 = 1:N
                    for l1 = 1:L
                        if n == n1 && l == l1
                            continue;
                        elseif n==n1 && A_opt(l1,n)==0 % case 2
                            A_new = A_opt; A_new(l,n) = 0; A_new(l1,n) = 1;
                            [cost_new, D, W] = optimize_D_matching(L, N, K, H_herm, G, G_diag, A_new, gamma_min, sigma2);
                            if cost_new < cost(end)
                                A_opt = A_new; D_opt = D; W_opt = W; cost=[cost, cost_new];
                            end
                        elseif A_opt(l1,n1) == 1 % case 3
                            A_new = A_opt; 
                            A_new(l,n) = 0; A_new(l1,n1) = 0; A_new(l,n1) = 1; A_new(l1,n) = 1; 
                            [cost_new, D, W] = optimize_D_matching(L, N, K, H_herm, G, G_diag, A_new, gamma_min, sigma2);
                            if cost_new < cost(end)
                                A_opt = A_new; D_opt = D; W_opt = W; cost=[cost, cost_new];
                            end
                        end 
                    end
                end
                % case 1
                A_new = A_opt; A_new(l,n) = 0;
                [cost_new, D, W] = optimize_D_matching(L, N, K, H_herm, G, G_diag, A_new, gamma_min, sigma2);
                if cost_new < cost(end)
                    A_opt = A_new; D_opt = D; W_opt = W; cost=[cost, cost_new];
                end
            end
        end
    end
    if abs(cost(end) - cost_t) <= 1e-4
        break;
    end
end
end
