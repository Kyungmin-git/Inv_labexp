load('Spect_store_300_figure.mat','spect_store_300');
load('Spect_store_50_figure.mat','f_h');
load('lsponge_noise.mat','spect_store_lsponge_noise');
load('Spect_store_dspongy.mat','spect_store_dspongy');

figure(1)
plot(f_h(1:1000),spect_store_300(2,1:1000),'LineWidth',1.5)
hold on
plot(f_h(1:1000),spect_store_dspongy(1:1000,5),'LineWidth',2,'Color',[0, 0.25, 0])
plot(f_h(1:1000),spect_store_300(1,1:1000),'LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980])


hold off
xlim([0 50])
ylim([0 120])
xlabel('Frequency (Hz)')
ylabel('Pressure (Pa*s)')
legend('Light sponge','Long Dense sponge','9-Capillary')
set(gca,'FontSize',18)
title('300 ml/s gas flow (Different porous media)')
