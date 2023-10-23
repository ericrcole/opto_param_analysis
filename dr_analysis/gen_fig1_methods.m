%generate figures of standard parameter spaces

addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/Framework'));
addpath('/Users/ERCOLE/Documents/Research/septum_review');
addpath('/Users/ERCOLE/Documents/Research/Repos/param_analysis');

fs = 500;
dur = 2;

tt1 = 0:1/fs:dur;
tt = 1/fs:1/fs:dur;
f = 1/fs:1/10:dur;

std = 20*pulstran(tt,f,.005,fs);

fig_position = [100 100 1200 800];
fsize = 18;

f1 = figure('Position',fig_position);
subplot(2,3,2)
plot(tt1,opto_generate_standard(fs,7,dur,20,.01,1),'k')
title('Standard pulse train')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
set(gca,'fontsize',fsize)
f1.Position = fig_position;
%saveas(gcf,'stim_standard.png')
% saveas(gcf,'stim_standard.fig')

subplot(2,3,3)
plot(tt,opto_generate_nested_train(fs, 7,40, dur, 20, .08, .005, 1),'k')
title('Nested pulse train')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
set(gca,'fontsize',fsize)
%saveas(gcf,'stim_npt.png')
% saveas(gcf,'stim_npt.fig')

subplot(2,3,5)
plot(tt,opto_generate_sinusoid_stim(fs, dur, 1 , 20, 7),'k')
title('Sinusoid')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
set(gca,'fontsize',fsize)
%saveas(gcf,'stim_sine.fig')

subplot(2,3,6)
plot(tt,opto_generate_sinusoid_stim(fs, dur, 2 , [20 .5], [7 30]),'k')
title('Double sinusoid')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
set(gca,'fontsize',fsize)

subplot(2,3,4)
psn = opto_generate_poisson(fs, 15, dur, 20 , .01, 1,.002);
plot(tt1,psn(1:end-1),'k')
title('Poisson pulse train')
xlabel('Time (s)')
ylabel('Amplitude (mW/mm^2)')
ylim([0 22.5])
set(gca,'fontsize',fsize)

annotation('textbox',[.08 .875 .1 .1],'String','a','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 32);
annotation('textbox',[.365 .875 .1 .1],'String','b','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 32);
annotation('textbox',[.645 .875 .1 .1],'String','c','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 32);
annotation('textbox',[.08 .42 .1 .1],'String','d','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 32);
annotation('textbox',[.365 .42 .1 .1],'String','e','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 32);
annotation('textbox',[.645 .42 .1 .1],'String','f','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 32);

saveas(gcf,'param_analysis_fig1.png')
saveas(gcf,'param_analysis_fig1.fig')