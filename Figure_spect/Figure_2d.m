load('Spect_store_100_store.mat','spect_store_100');
load('Spect_store_50ml.mat','f_h');

load('lsponge_noise.mat','spect_store_lsponge_noise');
load('Spect_store_dspongy.mat','spect_store_dspongy');


figure(1)
plot(f_h(1:1000),spect_store_100(1:1000,4),'LineWidth',1.5)
hold on
plot(f_h(1:1000),spect_store_lsponge_noise(1:1000,2),'Color',[0,0.24,0.54],'LineStyle','--','LineWidth', ...
    1.5)
plot(f_h(1:1000),spect_store_100(1:1000,2),'LineWidth',1.5,'Color',[0.4660, 0.6740, 0.1880])
plot(f_h(1:1000),spect_store_dspongy(1:1000,3),'LineWidth',2,'Color',[0, 0.25, 0])
plot(f_h(1:1000),spect_store_100(1:1000,3),'LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980])
plot(f_h(1:1000),spect_store_100(1:1000,1),'LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250])

hold off
xlim([0 50])
ylim([0 40])
xlabel('Frequency (Hz)')
ylabel('Pressure (Pa*s)')
legend('Light sponge','Light sponge (Noise)','Dense sponge','Long Dense sponge','9-Capillary','1-Capillary')
set(gca,'FontSize',18)
title('100 ml/s gas flow (Different porous media)')
