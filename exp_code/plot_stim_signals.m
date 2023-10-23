%generate figures of standard parameter spaces

addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/Framework'));
addpath('/Users/ERCOLE/Documents/Research/septum_review');

save_plots = false;

fs = 500;
dur = 1;

tt1 = 0:1/fs:dur;
tt = 1/fs:1/fs:dur;
f = 1/fs:1/10:dur;

std = 20*pulstran(tt,f,.001,fs);

fig_position = [100 100 400 320];


f1 = figure;
plot(tt1,opto_generate_standard(fs,7,dur,20,.01,1),'k')
suptitle('Standard Pulse Train')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
f1.Position = fig_position;
%saveas(gcf,'stim_standard.png')
if save_plots
saveas(gcf,'stim_standard.fig')
end

f2 = figure;
plot(tt,opto_generate_nested_train(fs, 7,40, dur, 20, .08, .005, 1),'k')
suptitle('Nested Pulse Train')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
f2.Position = fig_position;
%saveas(gcf,'stim_npt.png')
if save_plots
saveas(gcf,'stim_npt.fig')
end

f3 = figure;
plot(tt,opto_generate_sinusoid_stim(fs, dur, 1 , 20, 7),'k')
suptitle('Sinusoid')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
f3.Position = fig_position;
%saveas(gcf,'stim_sine.png')
if save_plots
saveas(gcf,'stim_sine.fig')
end

f4 = figure;
plot(tt,opto_generate_sinusoid_stim(fs, dur, 2 , [20 .5], [7 30]),'k')
suptitle('Double Sinusoid')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
f4.Position = fig_position;
%saveas(gcf,'stim_sine2n.png')
if save_plots
saveas(gcf,'stim_sine2n.fig')
end

psn = opto_generate_poisson(fs, 15, dur, 20 , .01, 1,.002);
f5 = figure;
plot(tt1,psn(1:end-1),'k')
suptitle('Poisson Pulse Train')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
f5.Position = fig_position;
%saveas(gcf,'stim_poisson.png')
if save_plots
saveas(gcf,'stim_poisson.fig')
end


%%  generate final version for methods figure 1 of param paper
save_plots = true;
lims = [-1 22];

figure('Position',[200 200 1200 300])
subplot(1,5,1)
plot(tt1,opto_generate_standard(fs,7,dur,20,.01,1),'k')
title('Standard Pulse')
ylabel('Amplitude (mW/mm^2)')
set(gca,'FontSize',16)
ylim(lims)

subplot(1,5,2)
plot(tt,opto_generate_nested_train(fs, 7,40, dur, 20, .08, .005, 1),'k')
title('Nested Pulse')
set(gca,'FontSize',16)
ylim(lims)

subplot(1,5,3)
psn = opto_generate_poisson(fs, 15, dur, 20 , .01, 1,.002);
plot(linspace(0,dur,length(psn)),psn,'k')
title('Poisson Pulse')
xlabel('Time (s)')
set(gca,'FontSize',16)
ylim(lims)

subplot(1,5,4)
plot(tt,opto_generate_sinusoid_stim(fs, dur, 1 , 20, 7),'k')
title('Sinusoid')
set(gca,'FontSize',16)
ylim(lims)

subplot(1,5,5)
plot(tt,opto_generate_sinusoid_stim(fs, dur, 2 , [20 .5], [7 30]),'k')
title('Double Sinusoid')
ylim(lims)

set(gca,'FontSize',16)

if save_plots
    saveas(gcf,'fig1_params.png')
    saveas(gcf,'fig1_params.svg')
end
