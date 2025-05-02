clear all; clc;
load('datasets\result\PASS_SU_gamma_min_L_12_5.mat');
figure1 = figure('PaperType','<custom>','PaperSize',[22 16],'Color',[1 1 1]);
plot(gamma_min_dB_all, mean(val_BnB_all,1), 'Marker','square','Color',[0 0 1], ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'PASS (BnB-Optimal)');hold on;
plot(gamma_min_dB_all, mean(val_BnB_equal,1), 'Marker','x','Color',[0.635294139385223 0.0784313753247261 0.184313729405403], ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'PASS (BnB-Equal)');hold on;
plot(gamma_min_dB_all, mean(val_SM_all,1), 'Marker','>','Color',[0 0.498039215803146 0],...
     'MarkerSize',12,'LineWidth',2.5,'DisplayName', 'PASS (Matching)');hold on;
plot(gamma_min_dB_all, mean(val_MIMO_all,1), 'Marker','o','Color',[0.929411768913269 0.694117665290833 0.125490203499794], ...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'MIMO');hold on;
plot(gamma_min_dB_all, mean(val_HybridMIMO_all,1), 'Marker','diamond','Color',[0 0 0],...
    'MarkerSize',12,'LineWidth',2.5, 'DisplayName', 'Massive MIMO');hold on;

xlabel('Minimum SINR requirement of user (dB)');
ylabel('Transmit power (dBm)');
set(gca,'FontSize',14);
legend(gca,'show');
set(legend,'Location','southeast');
set(gca, 'FontName', 'Arial'); 
savefig('figures\fig_SU_gamma.fig');
print(gcf, 'figures\fig_SU_gamma.pdf', '-dpdf', '-painters');

load('datasets\result\PASS_SU_gamma_min_L_12.mat');
fig2 = figure('PaperType','<custom>','PaperSize',[22 16],'Color',[1 1 1]);
plot(1:length(BnB_GUB_equal{1,6}),ones(1,length(BnB_GUB_equal{1,6}))*BnB_GUB_all{1,6}(end), ...
    'LineWidth', 2, 'Marker', '*', 'color', 'k', 'DisplayName', 'Optimal value, BnB-Optimal'); hold on;
plot(BnB_GLB_equal{1,6}, 'LineStyle', ':', ...
    'LineWidth', 2, 'Marker', '>', 'color', 'b', 'DisplayName', 'GLB, BnB-Equal'); hold on;
plot(BnB_GUB_equal{1,6}, 'LineStyle', '--', ...
    'LineWidth', 2, 'Marker', '>', 'color', 'b', 'DisplayName', 'GUB, BnB-Equal'); hold on;
plot(1:length(BnB_GUB_equal{1,6}),ones(1,length(BnB_GUB_equal{1,6}))*BnB_GUB_equal{1,6}(end), ...
    'LineStyle', '-', 'LineWidth', 2, 'Marker', 'o', 'color', 'b', 'DisplayName', 'Optimal value, BnB-Equal');
xlabel('Number of BnB iterations');
ylabel('Transmit power (dBm)');
set(gca,'FontSize',14);
legend(gca,'show');
set(legend,'interpreter','tex');
savefig('figures\fig_SU_BnBconv.fig');
print(gcf, 'figures\fig_SU_BnBconv.pdf', '-dpdf', '-painters');
