clear all; clc;
addpath(genpath(pwd));
datasetPath = strcat(pwd, '\datasets\');

S_x = 15;
S_y = 15;
K = 1;
N = 2;
freq = 15*10e9;
lightspeed = 3*10e8;
lambda = lightspeed/freq;
L = 12;
M = N*L;
PA_spacing = S_x/L;
regenerate = false;

noise_dBm = -80;
[H_conj_data, G, G_diag, sigma2, dist_PA_user, loc_U, n_samples] = cfgPinchingChannel(M,N,L,K,lambda,PA_spacing,freq,noise_dBm,lightspeed,S_x,S_y,regenerate);
[H_conj_data_MIMO, dist_MIMO_user] = cfgMIMOChannel(N, K, lambda, freq, noise_dBm, lightspeed, loc_U);
[H_conj_data_HybridMIMO, dist_HybridMIMO_user] = cfgHybridMIMOChannel(M, N, L, K, lambda, freq, noise_dBm, lightspeed, loc_U);

val_SM_all = []; val_BnB_all = []; val_BnB_equal = []; val_MIMO_all = []; val_HybridMIMO_all = [];
SM_conv_all = {}; BnB_GLB_all = {}; BnB_GUB_all = {}; BnB_GLB_equal = {}; BnB_GUB_equal = {}; BnB_GLB_equal_all = {}; BnB_GUB_equal_all = {};
gamma_min_dB_all = 20:1:30;
subconnected = true;
for idx = 1:10
    H_conj = H_conj_data(:,:,idx); % [M, K]
    H_herm = H_conj.'; % [K, M]
    H_herm_MIMO = H_conj_data_MIMO(:,:,idx).';
    H_herm_hybridMIMO = H_conj_data_HybridMIMO(:,:,idx).';
    dist_PU = dist_PA_user(:,:,idx); % [M, K]
    disp(['idx:',num2str(idx)]);
    parfor g = 1:length(gamma_min_dB_all)
        gamma_min_dB = gamma_min_dB_all(g);
        gamma_min = 10^(gamma_min_dB/10);

        [A_SM, D_SM, val_SM] =  Matching_SU(K, M, N, L, H_herm, G, G_diag, gamma_min, sigma2);
        val_SM_all(idx, g) = val_SM(end); SM_conv_all{idx,g} = val_SM;

        disp([num2str(g), 'SM completed']);

        [D_MIMO, val_MIMO] =  MIMO_SU(K, N, H_herm_MIMO, gamma_min, sigma2);
        val_MIMO_all(idx, g) = val_MIMO;

        disp([num2str(g), 'DMIMO completed']);

        [A_HybridMIMO, D_HybridMIMO, val_HybridMIMO] =  HybridMIMO_Penalty_SU(K, M, N, L, H_herm_hybridMIMO, gamma_min, sigma2, subconnected);
        val_HybridMIMO_all(idx, g) = val_HybridMIMO;

        disp([num2str(g), 'Hybrid completed']);
        
        val_BnB_best = +inf; val_BnB_equal_best = +inf;
        for l1 = 1:L
            for l2 = 1:L
                L_sel = [l1, l2];
                [A_BnB, D_BnB, val_BnB, GLB_hist, GUB_hist] =  BnB_SU(K, M, N, L, H_herm, G, G_diag, gamma_min, sigma2, L_sel);
                if val_BnB < val_BnB_best 
                    val_BnB_best = val_BnB;
                    val_BnB_all(idx, g) = val_BnB; 
                    BnB_GLB_all{idx, g} = GLB_hist; 
                    BnB_GUB_all{idx, g} = GUB_hist;
                end
                if l1 == l2
                    BnB_GLB_equal_all{idx, g, l1} = GLB_hist;
                    BnB_GUB_equal_all{idx, g, l1} = GUB_hist;
                end
                if l1==l2 && val_BnB < val_BnB_equal_best
                    val_BnB_equal_best = val_BnB;
                    val_BnB_equal(idx, g) = val_BnB;
                    BnB_GLB_equal{idx, g} = GLB_hist;
                    BnB_GUB_equal{idx, g} = GUB_hist;
                end
            end
        end
        disp([num2str(g), 'BnB completed']);        
    end
    save([datasetPath,'\result\PASS_SU_gamma_min_L_12.mat']);
end
