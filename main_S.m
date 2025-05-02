clear all; clc;
addpath(genpath(pwd));
datasetPath = strcat(pwd, '\datasets\');

freq = 15*10e9;
lightspeed = 3*10e8;
lambda = lightspeed/freq;
L_sel = -1;
Pmax_dBm = 10;
Pmax = dbm_to_watts(Pmax_dBm);

gamma_min_dB = 20;
gamma_min = 10^(gamma_min_dB/10);


regenerate = false;
subconnected = true;

S_y = 10;
S_x_all = 5:5:30;

for L = [10, 24]
    for K = [2, 4]
        N = K;
        M = N*L;
        val_SM_all = []; val_BnB_all = []; val_MIMO_all = []; val_HybridMIMO_all = []; 
        sol_SM_all = {};
        SM_conv_all = {}; BnB_GLB_all = {}; BnB_GUB_all = {}; 
        n_samples = 100;
        for idx = 1:20
            disp(['idx:',num2str(idx)]);
            parfor (s = 1:length(S_x_all), 6)
                S_x = S_x_all(s);
                PA_spacing = S_x/L;
                noise_dBm = -80;
                [H_conj_data, G, G_diag, sigma2, dist_PA_user, loc_U, n_samples] = cfgPinchingChannel(M, N,L,K,lambda,PA_spacing,freq,noise_dBm,lightspeed,S_x,S_y,regenerate);
                [H_conj_data_MIMO, dist_MIMO_user] = cfgMIMOChannel(N, K, lambda, freq, noise_dBm, lightspeed, loc_U);
                [H_conj_data_HybridMIMO, dist_HybridMIMO_user] = cfgHybridMIMOChannel(M, N, L, K, lambda, freq, noise_dBm, lightspeed, loc_U);
                H_conj = H_conj_data(:,:,idx); % [M, K]
                H_herm = H_conj.'; % [K, M]
                H_herm_MIMO = H_conj_data_MIMO(:,:,idx).';
                H_herm_hybridMIMO = H_conj_data_HybridMIMO(:,:,idx).';
                        
                [A_HybridMIMO, D_HybridMIMO, val_HybridMIMO] =  HybridMIMO_Penalty(K, M, N, L, H_herm_hybridMIMO, gamma_min, sigma2, subconnected);
                val_HybridMIMO_all(idx, s) = val_HybridMIMO;
        
                [A_SM, D_SM, val_SM] =  Matching(K, M, N, L, H_herm, G, G_diag, gamma_min, sigma2);
                val_SM_all(idx, s) = val_SM(end); SM_conv_all{idx,s} = val_SM;
                sol_SM_all{idx, s} =  struct('A', A_SM, 'D', D_SM); 
                
                [D_MIMO, val_MIMO] =  MIMO(K, N, H_herm_MIMO, gamma_min, sigma2);
                val_MIMO_all(idx, s) = val_MIMO;
            end
            save([datasetPath,'\result\PASS_S_K_',num2str(K),'_L_',num2str(L),'_withSol.mat']);
        end
    end
end