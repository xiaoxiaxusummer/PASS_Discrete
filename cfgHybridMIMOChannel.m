function [H_conj_data, dist_MIMO_user] = cfgHybridMIMOChannel(M,N,L,K,lambda,freq,noise_dBm,lightspeed,loc_U)
    h_MIMO = 5;

    antenna_spacing = lambda/2;
    
    kappa = 2*pi/lambda;
    noise = dbm_to_watts(noise_dBm);
    sigma2 = noise/noise; % 归一化的噪声功率
    
    beta = lightspeed/4/pi/freq;
    beta = beta/noise; % 归一化的信道增益
    
    
    loc_MIMO_antennas = zeros(L, N, 3); loc_MIMO_antennas(:,:,3) = h_MIMO;
    for n = 1:N 
        loc_MIMO_antennas(:,n,1) = 0:antenna_spacing:(L-1)*antenna_spacing;
        loc_MIMO_antennas(:,n,2) = (n-1)*antenna_spacing;
    end
    
    n_samples = 100;
    
    dist_MIMO_user = sqrt((reshape(loc_MIMO_antennas(:,:,1),[L,N,1,1]).*ones(L,N,K,n_samples) - reshape(loc_U(:,:,1),[1,1,K,n_samples]).*ones(L,N,K,n_samples)).^2 + ...
                    (reshape(loc_MIMO_antennas(:,:,2),[L,N,1,1]).*ones(L,N,K,n_samples) - reshape(loc_U(:,:,2),[1,1,K, n_samples]).*ones(L,N,K,n_samples)).^2 + ... 
                    (h_MIMO)^2); % [L, N, K, n_samples]
    H_conj_data = sqrt(beta) ./dist_MIMO_user .* exp(-1j*kappa*dist_MIMO_user); % [L, N, K, n_samples]
    % Matlab优先逐列拼接
    H_conj_data = reshape(H_conj_data,[M,K,n_samples]); % [M, K, n_samples]
end