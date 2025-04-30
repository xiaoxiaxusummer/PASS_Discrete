function [H_conj_data, dist_MIMO_user] = cfgMIMOChannel(N,K,lambda,freq,noise_dBm,lightspeed,loc_U)
    H = 5;

    antenna_spacing = lambda/2;
    
    kappa = 2*pi/lambda;

    noise = dbm_to_watts(noise_dBm);
    
    beta = lightspeed/4/pi/freq;
    beta = beta/noise; % 归一化的信道增益

    loc_MIMO_antennas = zeros(N, 3); loc_MIMO_antennas(:,3) = H;
    loc_MIMO_antennas(:,1) = 0:antenna_spacing:(N-1)*antenna_spacing;
    loc_MIMO_antennas(:,2) = (N-1)*antenna_spacing;
    
    n_samples = 100;
    dist_MIMO_user = sqrt((reshape(loc_MIMO_antennas(:,1),[N,1,1]).*ones(N,K,n_samples) - reshape(loc_U(:,:,1),[1,K,n_samples]).*ones(N,K,n_samples)).^2 + ...
                    (reshape(loc_MIMO_antennas(:,2),[N,1,1]).*ones(N,K,n_samples) - reshape(loc_U(:,:,2),[1,K, n_samples]).*ones(N,K,n_samples)).^2 + ... 
                    (H)^2); % [N, K, n_samples]
    H_conj_data = sqrt(beta) ./dist_MIMO_user .* exp(-1j*kappa*dist_MIMO_user); % [N, K, n_samples]
end