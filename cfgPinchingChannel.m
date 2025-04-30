function [H_conj_data, G, G_diag, sigma2, dist_PA_user, loc_U, n_samples] = cfgPinchingChannel(M,N,L,K,lambda,PA_spacing,freq,noise_dBm,lightspeed,S_x,S_y,regenerate)
h_PAA = 5;
n_eff = 1.4;

noise = dbm_to_watts(noise_dBm);
sigma2 = noise/noise; % 归一化的噪声功率

kappa = 2*pi/lambda;
beta = lightspeed/4/pi/freq;
beta = beta/noise; % Normalized channel gain


loc_PA = zeros(L, N, 3); loc_PA(:,:,3) = h_PAA;
loc_W = zeros(N, 3); loc_W(:,3) = h_PAA;
if K==1
    for n = 1:N
        loc_PA(:,n,1) = PA_spacing:PA_spacing:L*PA_spacing; loc_PA(:,n,2) = n*S_y/2;
        loc_W(n,1:2) = [0, n*S_y/2];
    end
else
    for n = 1:N
        if n <= N/2
            loc_PA(:,n,1) = 0:PA_spacing:(L-1)*PA_spacing;
            loc_PA(:,n,2) = (n-1)*S_y/2;
            loc_W(n,1) = 0;
            loc_W(n,2) = (n-1)*S_y/2;
        else
            loc_PA(:,n,1) = (n-1)*S_x/2;
            loc_PA(:,n,2) = 0:PA_spacing:(L-1)*PA_spacing;
            loc_W(n,1) = (n-1)*S_x/2;
            loc_W(n,2) = 0;
        end
    end
end

loc_PA = reshape(loc_PA, [M,3]);

n_samples = 100;

if regenerate
    loc_U = zeros(K, n_samples, 3);
    loc_U(:,:,1) = rand(K,n_samples) * S_x/2 + S_x/2;
    loc_U(:,:,2) = rand(K,n_samples) * S_y/2 + S_y/2;
    save(['datasets\20250330_PASS_K_',num2str(K),'.mat']);
else
    loc_data = load(['datasets\20250330_PASS_K_',num2str(K),'.mat'],'loc_U','S_x','S_y');
    loc_U = loc_data.loc_U; S_x0 = loc_data.S_x; S_y0 = loc_data.S_y;
    loc_U(:,:,1) = loc_U(:,:,1) / S_x0 * S_x;
    loc_U(:,:,2) = loc_U(:,:,2) / S_y0 * S_y;
end

dist_PA_user = sqrt((reshape(loc_PA(:,1),[M,1,1]).*ones(M,K,n_samples) - reshape(loc_U(:,:,1),[1,K,n_samples]).*ones(M,K,n_samples)).^2 + ...
    (reshape(loc_PA(:,2),[M,1,1]).*ones(M,K,n_samples) - reshape(loc_U(:,:,2),[1,K,n_samples]).*ones(M,K,n_samples)).^2 + ...
    (h_PAA)^2); % [L, N, K, n_samples]

dist_PA_waveguide = sqrt((reshape(loc_PA(:,1),[L,N]) - loc_W(:,1).'.*ones(L,N)).^2 + (reshape(loc_PA(:,2),[L,N])  - loc_W(:,2).'.*ones(L,N)).^2); % [L,N]

G_ori = exp(-1j*kappa*n_eff*dist_PA_waveguide); % [L, N]
G = zeros(M, N);
for n = 1:N
    G((n-1)*L+1:n*L, n) = G_ori(:, n);
end
G_diag = diag(reshape(G_ori, [M,1]));

H_conj_data = sqrt(beta) ./dist_PA_user .* exp(-1j*kappa*dist_PA_user); % [M, K, n_samples]
end