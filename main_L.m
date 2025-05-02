clear all; clc;
addpath(genpath(pwd));
datasetPath = strcat(pwd, '\datasets\');

S_x = 10;
S_y = 10;


freq = 15*10e9;
lightspeed = 3*10e8;
lambda = lightspeed/freq;
L_sel = -1;


Pmax_dBm = 10;
Pmax = dbm_to_watts(Pmax_dBm);

gamma_min_dB = 20;
gamma_min = 10^(gamma_min_dB/10);

regenerate = false;
L_all = 4:2:20;


subconnected = true;
n_samples = 100;

for K = [2, 4]
    for S_x = [10]
        N = K;
        val_SM_all = []; val_BnB_all = []; val_MIMO_all = []; val_HybridMIMO_all = [];
        SM_conv_all = {}; BnB_GLB_all = {}; BnB_GUB_all = {}; 
        sol_SM_all = {};
        for idx = 1:20
            disp(['idx:',num2str(idx)]);
            parfor (l = 1:length(L_all),11)
                L = L_all(l);
                M = N*L;
                PA_spacing = S_x/L;
                noise_dBm = -80;
                [H_conj_data, G, G_diag, sigma2, dist_PA_user, loc_U, n_samples] = cfgPinchingChannel(M, N,L,K,lambda,PA_spacing,freq,noise_dBm,lightspeed,S_x,S_y,regenerate);
                [H_conj_data_MIMO, dist_MIMO_user] = cfgMIMOChannel(N, K, lambda, freq, noise_dBm, lightspeed, loc_U);
                [H_conj_data_HybridMIMO, dist_HybridMIMO_user] = cfgHybridMIMOChannel(M, N, L, K, lambda, freq, noise_dBm, lightspeed, loc_U);
                H_conj = H_conj_data(:,:,idx); % [M, K]
                H_herm = H_conj.'; % [K, M]
                H_herm_MIMO = H_conj_data_MIMO(:,:,idx).';
                H_herm_hybridMIMO = H_conj_data_HybridMIMO(:,:,idx).';
                dist_PU = dist_PA_user(:,:,idx); % [L, N, K]
                [A_HybridMIMO, D_HybridMIMO, val_HybridMIMO] =  HybridMIMO_Penalty(K, M, N, L, H_herm_hybridMIMO, gamma_min, sigma2, subconnected);
                val_HybridMIMO_all(idx, l) = val_HybridMIMO;
        
                [A_SM, D_SM, val_SM] =  Matching(K, M, N, L, H_herm, G, G_diag, gamma_min, sigma2);
                val_SM_all(idx, l) = val_SM(end); SM_conv_all{idx,l} = val_SM;
                sol_SM_all{idx, l} = struct("A", A_SM, "D", D_SM); 
        
                [D_MIMO, val_MIMO] =  MIMO(K, N, H_herm_MIMO, gamma_min, sigma2);
                val_MIMO_all(idx, l) = val_MIMO;
            end
            save([datasetPath,'\result\PASS_L_K_',num2str(K),'_Sx_',num2str(S_x),'.mat']);
        end
    end
end