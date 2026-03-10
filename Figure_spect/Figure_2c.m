load('Spect_store_50ml_new.mat','f_h','spect_store');

spect_store_50=spect_store;
load('lsponge_noise.mat','spect_store_lsponge_noise');
load('Spect_store_dspongy.mat','Spect_store');
spect_store_dspongy50=Spect_store;

figure(1)
plot(f_h(1:1000),spect_store_50(1,1:1000),'LineWidth',1.5)
hold on
plot(f_h(1:1000),spect_store_lsponge_noise(1:1000,1),'Color',[0,0.24,0.54],'LineStyle','--','LineWidth', ...
    1.5)
plot(f_h(1:1000),spect_store_50(2,1:1000),'LineWidth',1.5,'Color',[0.4660, 0.6740, 0.1880])
plot(f_h(1:1000),spect_store_dspongy50(1:1000,2),'LineWidth',2,'Color',[0, 0.25, 0])
plot(f_h(1:1000),spect_store_50(3,1:1000),'LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980])
plot(f_h(1:1000),spect_store_50(4,1:1000),'LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250])

hold off
xlim([0 50])
ylim([0 25])
xlabel('Frequency (Hz)')
ylabel('Pressure (Pa*s)')
legend('Light sponge','Light sponge (Noise)','Dense sponge','Long Dense sponge','9-Capillary','1-Capillary')
set(gca,'FontSize',18)
title('50 ml/s gas flow (Different porous media)')
