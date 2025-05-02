clear all; clc;
addpath(genpath(pwd));
datasetPath = strcat(pwd, '\datasets\');

S_x = 10;
S_y = 10;

K = 2;
N = K;
freq = 15*10e9;
lightspeed = 3*10e8;
lambda = lightspeed/freq;
L = 6;
M = N*L;
PA_spacing = S_x/L;


P0_dBm = -10;
P0 = dbm_to_watts(P0_dBm);

regenerate = false;

val_SM_all = []; val_BnB_all = []; val_MIMO_all = []; val_HybridMIMO_all = [];
SM_conv_all = {}; BnB_GLB_all = {}; BnB_GUB_all = {}; 


noise_dBm = -80;
[H_conj_data, G, G_diag, sigma2, dist_PA_user, loc_U, n_samples] = cfgPinchingChannel(M, N,L,K,lambda,PA_spacing,freq,noise_dBm,lightspeed,S_x,S_y,regenerate);
[H_conj_data_MIMO, dist_MIMO_user] = cfgMIMOChannel(N, K, lambda, freq, noise_dBm, lightspeed, loc_U);
[H_conj_data_HybridMIMO, dist_HybridMIMO_user] = cfgHybridMIMOChannel(M, N, L, K, lambda, freq, noise_dBm, lightspeed, loc_U);


gamma_min_dB_all = 20:1:30;
subconnected = true;
for idx = 1:20
    H_conj = H_conj_data(:,:,idx); % [M, K]
    H_herm = H_conj.'; % [K, M]
    H_herm_MIMO = H_conj_data_MIMO(:,:,idx).';
    H_herm_hybridMIMO = H_conj_data_HybridMIMO(:,:,idx).';
    dist_PU = dist_PA_user(:,:,idx); % [L, N, K]
    disp(['idx:',num2str(idx)]);
    parfor g = 1:length(gamma_min_dB_all)
        gamma_min_dB = gamma_min_dB_all(g);
        gamma_min = 10^(gamma_min_dB/10);
        
        [A_BnB, D_BnB, val_BnB, GLB_hist, GUB_hist] =  BnB(K, M, N, L, H_herm, G, G_diag, P0, gamma_min, sigma2);
        val_BnB_all(idx, g) = val_BnB; 
        BnB_GLB_all{idx, g} = GLB_hist; 
        BnB_GUB_all{idx, g} = GUB_hist;
        disp([num2str(g), 'BnB completed']);

        [A_SM, D_SM, val_SM] =  Matching(K, M, N, L, H_herm, G, G_diag, gamma_min, sigma2);
        val_SM_all(idx, g) = val_SM(end); SM_conv_all{idx,g} = val_SM;

        disp([num2str(g), 'Matching completed']);
        
        [D_MIMO, val_MIMO] =  MIMO(K, N, H_herm_MIMO, gamma_min, sigma2);
        val_MIMO_all(idx, g) = val_MIMO;

        disp([num2str(g), 'MIMO completed']);
        
        [A_HybridMIMO, D_HybridMIMO, val_HybridMIMO] =  HybridMIMO_Penalty(K, M, N, L, H_herm_hybridMIMO, gamma_min, sigma2, subconnected);
        val_HybridMIMO_all(idx, g) = val_HybridMIMO;

        disp([num2str(g), 'Hybrid MIMO completed']);
    end
    save([datasetPath,'\result\PASS_gamma_min.mat']);
end



function parsave(filePath, idx, gamma_min_dB, val_SM, SM_conv, val_MIMO, val_HybridMIMO, val_BnB, BnB_GLB, BnB_GUB, ...
    noise_dBm, L, K, M, N, S_x, S_y, PA_spacing, freq, lambda, subconnected)
    save(filePath);
end