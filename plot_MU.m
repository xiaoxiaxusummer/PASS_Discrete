clear all; clc;
%% P - Gamma_min
load('datasets\result\PASS_gamma_min.mat');
figure1 = figure('PaperType','<custom>','PaperSize',[22 16],'Color',[1 1 1]);
axes = axes('Parent',figure1,'Position',[0.13 0.119205298013245 0.775 0.805794701986755]);
hold(axes,'on');
plot(gamma_min_dB_all, mean(val_BnB_all,1), 'Marker','square','Color',[0 0 1], ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'PASS (BnB)');hold on;
plot(gamma_min_dB_all, mean(val_SM_all,1), 'Marker','>','Color',[0 0.498039215803146 0],'LineStyle','--',...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Matching)');hold on;
plot(gamma_min_dB_all, mean(val_MIMO_all,1), 'Marker','o','Color',[0.929411768913269 0.694117665290833 0.125490203499794], ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'MIMO');hold on;
plot(gamma_min_dB_all, mean(val_HybridMIMO_all,1), 'Marker','diamond','Color',[0 0 0],...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'Massive MIMO');hold on;

xlabel('Minimum rate requirement of user (dB)');
ylabel('Transmit power (dBm)');
box(axes,'on');
hold(axes,'off');
set(axes,'FontSize',14);
legend(gca,'show','interpreter','latex');
savefig('figures\fig_MU_gamma.fig');
print(gcf, 'figures\fig_MU_gamma.pdf', '-dpdf', '-painters');

%% Convergence of BnB
fig2 = figure('PaperType','<custom>','PaperSize',[22 16],'Color',[1 1 1]);
plot(BnB_GLB_all{4,1},'DisplayName', 'GLB'); hold on;
plot(BnB_GUB_all{4,1},'DisplayName', 'GUB'); hold on;
plot(0:100:1500,ones(1,length(0:100:1500))*BnB_GUB_all{4,1}(end),'DisplayName', 'Optimum');
xlabel('Number of BnB iterations');
ylabel('Transmit power (dBm)');
set(gca,'FontSize',14);
legend(gca,'show','interpreter','tex');
savefig('figures\fig_MU_BnBconv.fig');
print(gcf, 'figures\fig_MU_BnBconv.pdf', '-dpdf', '-painters');

%% Convergence of matching
fig2 = figure('PaperType','<custom>','PaperSize',[22 16],'Color',[1 1 1]);
plot(SM_conv_all{4,1},'DisplayName', 'Many-to-many Matching'); hold on;
plot(ones(1,length(SM_conv_all{4,1}))*BnB_GUB_all{4,1}(end),'DisplayName', 'Optimum');
xlabel('Number of BnB iterations');
ylabel('Transmit power (dBm)');
set(gca,'FontSize',14);
legend(gca,'show','interpreter','tex');
savefig('figures\fig_MU_SMconv.fig');
print(gcf, 'figures\fig_MU_SMconv.pdf', '-dpdf', '-painters');


fig3 = figure('PaperType','<custom>','PaperSize',[22 16],'Color',[1 1 1]);
load('datasets\result\PASS_L_K_2_Sx_10.mat');
plot(L_all, mean(val_SM_all,1), 'Marker','>','Color',[0 0.498039215803146 0],'LineStyle','--',...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Matching), $K=2$');hold on;
plot(L_all, mean(val_MIMO_all,1), 'Marker','o','Color',[0.929411768913269 0.694117665290833 0.125490203499794], 'LineStyle','--', ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'MIMO, $K=2$');hold on;
plot(L_all, mean(val_HybridMIMO_all,1), 'Marker','diamond','Color',[0 0 0], 'LineStyle','--', ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'Massive MIMO, $K=2$');hold on;
load('datasets\result\20250404_PASS_L_K_4_Sx_10_withSol.mat');
plot(L_all, mean(val_SM_all,1), 'Marker','>','Color',[0 0.498039215803146 0],...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Matching), $K=4$');hold on;
plot(L_all, mean(val_MIMO_all,1), 'Marker','o','Color',[0.929411768913269 0.694117665290833 0.125490203499794], ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'MIMO, $K=4$');hold on;
plot(L_all, mean(val_HybridMIMO_all,1), 'Marker','diamond','Color',[0 0 0],...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'Massive MIMO, $K=4$');hold on;
legend(gca,'show','interpreter','latex');
set(gca, 'FontName', 'Arial'); 
savefig('figures\fig_MU_L.fig');
print(gcf, 'figures\fig_MU_L.pdf', '-dpdf', '-painters');

fig4 = figure('PaperType','<custom>','PaperSize',[22 16],'Color',[1 1 1]);
load('datasets\result\continuousPASS_S_K_2_L_10.mat');
plot(S_x_all, mean(val_continuous_all,1), 'Marker','o','Color',[0 0.498039215803146 0],'LineStyle','--',...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Continuous), $K=2$');hold on;
load('datasets\result\PASS_S_K_2_L_10.mat');
plot(S_x_all, mean(val_SM_all,1), 'Marker','o','Color',[0 0.498039215803146 0],'LineStyle','--',...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Discrete), $K=2$');hold on;
plot(S_x_all, mean(val_HybridMIMO_all,1), 'Marker','diamond','Color',[0 0 0],...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'Massive MIMO, $K=2$');hold on;
load('datasets\result\continuousPASS_S_K_4_L_10.mat');
plot(S_x_all, mean(val_continuous_all,1), 'Marker','o','Color',[0 0.498039215803146 0],'LineStyle','--',...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Continuous), $K=4$');hold on;
load('datasets\result\PASS_S_K_4_L_10.mat');
plot(S_x_all, mean(val_SM_all,1), 'Marker','o','Color',[0 0.498039215803146 0],'LineStyle','--',...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Discrete), $K=4$');hold on;
plot(S_x_all, mean(val_HybridMIMO_all,1), 'Marker','diamond','Color',[0 0 0],...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'Massive MIMO, $K=4$');hold on;
xlabel("Spatial range (meter)");
ylabel("Transmit power (dBm)");
set(gca,'FontSize',14);
legend(gca,'show','interpreter','latex');
set(gca, 'FontName', 'Arial');
savefig('figures\fig_MU_S.fig');
print(gcf, 'figures\fig_MU_S.pdf', '-dpdf', '-painters');

