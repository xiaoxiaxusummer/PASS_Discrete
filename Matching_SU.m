function [A_opt, d_opt, cost] =  Matching_SU(K, M, N, L, H_herm, G, G_diag, gamma_min, sigma2)
gain_mean = reshape(mean(abs(H_herm).^2,1).',[L,N]);
A0 = zeros(L,N);
% g_max = -inf;
% for n = 1:N
%     [g_max_n,ind]=max(gain_mean(:,n));
%     if g_max_n > g_max
%         g_max = g_max_n; A0 = zeros(L,N);  A0(ind,n) = 1;
%     end
% end
% for n = 1:N
%     [~,ind]=max(gain_mean(:,n));
%     A0(ind,n) = 1;
% end
A_opt = A0;
% [cost, d_opt] = optimize_d_matching_SU(L, N, A_opt, H_herm, G, gamma_min, sigma2);
cost = +inf;
cost = [cost]; iter = 0;
% Welfare-driven matching
while true
    cost_t = cost(end); iter = iter + 1;
    for n = 1:N
        for l= 1:L
            if A_opt(l,n) == 0
                % From (l,0)&(n,0) -> (l,n)&(0,0)
                A_new = A_opt; A_new(l,n) = 1;
                [cost_new, d] = optimize_d_matching_SU(L, N, A_new, H_herm, G, gamma_min, sigma2);
                if cost_new < cost
                    A_opt = A_new; d_opt = d; cost=[cost, cost_new];
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
                        elseif A_opt(l1,n1) == 1 % case 3
                            A_new = A_opt;
                            A_new(l,n) = 0; A_new(l1,n1) = 0; A_new(l,n1) = 1; A_new(l1,n) = 1;
                            [cost_new, d] = optimize_d_matching_SU(L, N, A_new, H_herm, G, gamma_min, sigma2);
                            if cost_new < cost(end)
                                A_opt = A_new; d_opt = d; cost=[cost, cost_new];
                            end
                        elseif n==n1 && A_opt(l1,n)==0 % case 2
                            A_new = A_opt; A_new(l,n) = 0; A_new(l1,n) = 1;
                            [cost_new, d] = optimize_d_matching_SU(L, N, A_new, H_herm, G, gamma_min, sigma2);
                            if cost_new < cost(end)
                                A_opt = A_new; d_opt = d; cost=[cost, cost_new];
                            end

                        end
                    end
                end
                % case 1
                A_new = A_opt; A_new(l,n) = 0;
                [cost_new, d] = optimize_d_matching_SU(L, N, A_new, H_herm, G, gamma_min, sigma2);
                if cost_new < cost(end)
                    A_opt = A_new; d_opt = d; cost=[cost, cost_new];
                end
            end
        end
    end
    if abs(cost(end) - cost_t) <= 1e-4
        break;
    end
end
end

function [obj, d_opt] = optimize_d_matching_SU(L, N, A, H_herm, G, gamma_min, sigma2)
M = N*L;
L_sel = sum(A,1);
effective_gain = abs(H_herm * diag(A(:)) * G).^2;
activated_gain = 0;
for n = 1:N
    if L_sel(n) > 0
        activated_gain = activated_gain +  effective_gain(n)./ L_sel(n);
    end
end
% P_opt = gamma_min * sigma2 / sum(abs(H_herm * diag(A(:)) * G).^2 ./ L_sel);
P_opt = gamma_min * sigma2 / activated_gain;
d_opt = sqrt(P_opt) * (H_herm * diag(A(:)) * G)' / sqrt(sum(abs(H_herm * diag(A(:)) * G).^2));
obj = watts_to_dbm(P_opt);
end

