clear all; clc;
addpath(genpath(pwd));
currentFolder = fileparts(mfilename('fullpath'));
S_x_all = 5:5:30; S_y = 10;
L_all = [10];
for K = [4, 2]
    for L = L_all
        val_continuous_all = []; sol_continuous_all = {};
        for idx = 1:20
            parfor s = 1:length(S_x_all)
                N = K;
                M = N*L;
                S_x = S_x_all(s);
                cfg = load(['datasets\20250330_PASS_K_',num2str(K),'.mat']);
                loc_U = cfg.loc_U;
                if ~(cfg.S_x == S_x) || ~(cfg.S_y == S_y)
                    loc_U(:,:,1) = loc_U(:,:,1)/cfg.S_x * S_x;
                    loc_U(:,:,2) = loc_U(:,:,2)/cfg.S_y * S_y;
                end
                A = loc_U(:,:,1);
                cfg.S_x = S_x;
                cfg.L = L; cfg.N = N; cfg.K = K; cfg.M = M;
                cfg.sigma2 = cfg.sigma2;
                cfg.gamma_min_dB = 20; cfg.gamma_min = 10^(cfg.gamma_min_dB/10);
                cfg.loc_PA = zeros(M,3);
                cfg.loc_PA(:,1) = reshape((0:1:L-1).'*cfg.lambda/2 .* ones(L,N), [M,1]);
                cfg.loc_PA(:,2) = reshape((0:1:N-1)* (S_y/N) .* ones(L,N), [M,1]);
                cfg.loc_PA(:,3) = ones(M,1)*cfg.h_PAA;
                cfg.loc_W = zeros(N,3);
                cfg.loc_W(:,2) = (0:1:N-1)* (S_y/N);
                cfg.loc_W(:,3) = ones(N,1)*cfg.h_PAA;
                loc_U_data = loc_U; A_data = A;
                loc_U = reshape(loc_U_data(:,idx,:),[K,3]);
                A = A_data(:,idx).';
                phi = sqrt( (cfg.loc_PA(:,2) - loc_U(:,2).').^2 + cfg.h_PAA^2); % [M,N]
                discrete_data = load(['datasets\result\PASS_S_K_',num2str(K),'_L_',num2str(L),'.mat']);
                discrete_sol = discrete_data.sol_SM_all{idx,s};
                discrete_spacing = S_x/L;
                discrete_A = discrete_sol.A; discrete_D = discrete_sol.D; discrete_W = discrete_D.*sqrt(sum(discrete_A,1).');
                discrete_x = (0:discrete_spacing:discrete_spacing*(L-1)).' .* ones(L,N);
                [obj,sol,power_log,residuals] = run_continuous_search(M, N, L, K, loc_U, A, phi, cfg, discrete_A, discrete_x, discrete_W);
                val_continuous_all(idx, s) = obj; sol_continuous_all{idx, s} = sol;
            end
        end
        save(['datasets\result\continuousPASS_S_K_',num2str(K),'_L_', num2str(L),'.mat']);
    end
end

