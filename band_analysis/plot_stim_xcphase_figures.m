%param analysis - generating figures of autocorrelation for different parameter spaces 
%for thomas, along with violin plot of phase stuff

%cd('/Users/ERCOLE/Downloads/OneDrive_1_3-11-2022')
cd('/Volumes/OPTO_PROC/OptoStim')

subject_IDs = {'STV003','STV008'};

addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/param_spaces_analysis/Violinplot-Matlab'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/param_analysis'));

theta_ac_std = [];
gamma_ac_std =  [];

theta_ac_sine =  [];
gamma_ac_sine =  [];

theta_ac_npt =  [];
gamma_ac_npt =  [];

theta_ac_psn =  [];
gamma_ac_psn =  [];

theta_ac_rec =  [];
gamma_ac_rec =  [];

theta_phi_std =  [];
gamma_phi_std =  [];
theta_phi_sine =  [];
gamma_phi_sine =  [];
theta_phi_npt =  [];
gamma_phi_npt =  [];
theta_phi_psn =  [];
gamma_phi_psn =  [];
theta_phi_rec =  [];
gamma_phi_rec =  [];

for k = 1:length(subject_IDs)
    
    subject_ID = subject_IDs{k};
    if strcmp(subject_ID,'STV008')
        amp_thr = 75;
    else
        amp_thr = 30;
    end
        
    load(sprintf('%s_BandpowerPhase_Awake_Standard',subject_ID));
    bp_std = bandpower_data;
    load(sprintf('%s_BandpowerPhase_Awake_Sine',subject_ID));
    bp_sine = bandpower_data;
    load(sprintf('%s_BandpowerPhase_Awake_NPT',subject_ID));
    bp_npt = bandpower_data;
    load(sprintf('%s_BandpowerPhase_Awake_Poisson',subject_ID));
    bp_psn = bandpower_data;
    load(sprintf('%s_BandpowerPhase_Awake_Recording',subject_ID));
    bp_rec = bandpower_data;
    %filt_notch = 1;
    
    FSd = 2000;
    
    f_theta = 7;
    f_gamma = 35;
    
    
    fr_inds = bp_std.phi.fr_inds;
    theta_ind = find(fr_inds>=f_theta,1);
    gamma_ind = find(fr_inds>=f_gamma,1);
    
    
    % f_theta = 37;
    % f_gamma = 177;
    
    theta_inds_std = find(and(bp_std.param(:,2)==7,bp_std.param(:,1)>=amp_thr));
    theta_inds_sine = find(and(bp_sine.param(:,2)==7,bp_sine.param(:,1)>=amp_thr));
    %theta_inds_npt = find(and(bp_npt.param(:,2)==7,and(bp_npt.param(:,3)==50,bp_npt.param(:,1)>=amp_thr)));
    theta_inds_npt = find(and(and(bp_npt.param(:,2)>6,bp_npt.param(:,2)<8),bp_npt.param(:,1)>=amp_thr));
    theta_inds_psn = find(and(bp_psn.param(:,2)==7,bp_psn.param(:,1)>=amp_thr));
    
    gamma_inds_std = find(and(bp_std.param(:,2)==35,bp_std.param(:,1)>=amp_thr));
    gamma_inds_sine = find(and(bp_sine.param(:,2)==35,bp_sine.param(:,1)>=amp_thr));
    %gamma_inds_npt = find(and(bp_npt.param(:,2)==11,and(bp_npt.param(:,3)==35,bp_npt.param(:,1)>=amp_thr)));
    gamma_inds_npt = find(and(and(bp_npt.param(:,3)>34,bp_npt.param(:,3)<36),bp_npt.param(:,1)>=amp_thr));
    gamma_inds_psn = find(and(bp_psn.param(:,2)==35,bp_psn.param(:,1)>=amp_thr));
    
    if strcmp(subject_ID,'ARN088')
            %excluding trials from one experiment that had a really weird
            %~180 Hz noise harmonic that was corrupting the autocorr signal
        theta_inds_sine((theta_inds_sine<126)) = [];
        gamma_inds_sine((gamma_inds_sine<126)) = [];
        theta_inds_sine(and(theta_inds_sine>275,theta_inds_sine<663)) = [];
        gamma_inds_sine(and(gamma_inds_sine>275,gamma_inds_sine<663)) = [];
        
    end
    
    theta_ac_std = [theta_ac_std; cell2mat(cellfun(@transpose,bp_std.autocorr.stim(theta_inds_std),'UniformOutput',false))];
    gamma_ac_std = [gamma_ac_std; cell2mat(cellfun(@transpose,bp_std.autocorr.stim(gamma_inds_std),'UniformOutput',false))];
    
    theta_ac_sine = [theta_ac_sine; cell2mat(cellfun(@transpose,bp_sine.autocorr.stim(theta_inds_sine),'UniformOutput',false))];
    gamma_ac_sine = [gamma_ac_sine; cell2mat(cellfun(@transpose,bp_sine.autocorr.stim(gamma_inds_sine),'UniformOutput',false))];
    
    theta_ac_npt = [theta_ac_npt; cell2mat(cellfun(@transpose,bp_npt.autocorr.stim(theta_inds_npt),'UniformOutput',false))];
    gamma_ac_npt = [gamma_ac_npt; cell2mat(cellfun(@transpose,bp_npt.autocorr.stim(gamma_inds_npt),'UniformOutput',false))];
    
    theta_ac_psn = [theta_ac_psn; cell2mat(cellfun(@transpose,bp_psn.autocorr.stim(theta_inds_psn),'UniformOutput',false))];
    gamma_ac_psn = [gamma_ac_psn; cell2mat(cellfun(@transpose,bp_psn.autocorr.stim(gamma_inds_psn),'UniformOutput',false))];
    
    theta_ac_rec = [theta_ac_rec; cell2mat(cellfun(@transpose,bp_rec.autocorr.stim,'UniformOutput',false))];
    gamma_ac_rec = [gamma_ac_rec; cell2mat(cellfun(@transpose,bp_rec.autocorr.stim,'UniformOutput',false))];
    
    theta_phi_std = [theta_phi_std; cell2mat(cellfun(@(c)c(theta_ind),bp_std.phi.stim(theta_inds_std),'UniformOutput',false))];
    gamma_phi_std = [gamma_phi_std; cell2mat(cellfun(@(c)c(gamma_ind),bp_std.phi.stim(gamma_inds_std),'UniformOutput',false))];
    theta_phi_sine = [theta_phi_sine; cell2mat(cellfun(@(c)c(theta_ind),bp_sine.phi.stim(theta_inds_sine),'UniformOutput',false))];
    gamma_phi_sine = [gamma_phi_sine; cell2mat(cellfun(@(c)c(gamma_ind),bp_sine.phi.stim(gamma_inds_sine),'UniformOutput',false))];
    theta_phi_npt = [theta_phi_npt; cell2mat(cellfun(@(c)c(theta_ind),bp_npt.phi.stim(theta_inds_npt),'UniformOutput',false))];
    gamma_phi_npt = [gamma_phi_npt; cell2mat(cellfun(@(c)c(gamma_ind),bp_npt.phi.stim(gamma_inds_npt),'UniformOutput',false))];
    theta_phi_psn = [theta_phi_psn; cell2mat(cellfun(@(c)c(theta_ind),bp_psn.phi.stim(theta_inds_psn),'UniformOutput',false))];
    gamma_phi_psn = [gamma_phi_psn; cell2mat(cellfun(@(c)c(gamma_ind),bp_psn.phi.stim(gamma_inds_psn),'UniformOutput',false))];
    theta_phi_rec = [theta_phi_rec; cell2mat(cellfun(@(c)c(theta_ind),bp_rec.phi.stim,'UniformOutput',false))];
    gamma_phi_rec = [gamma_phi_rec; cell2mat(cellfun(@(c)c(gamma_ind),bp_rec.phi.stim,'UniformOutput',false))];
end

theta_phi_vplot = padcat(theta_phi_sine,theta_phi_std,theta_phi_psn,theta_phi_npt,theta_phi_rec);
gamma_phi_vplot = padcat(gamma_phi_sine,gamma_phi_std,gamma_phi_psn,gamma_phi_npt,gamma_phi_rec);


%cellfun(@(c,idx)c(idx,:),C,

%mean, standard deviation, and standard error across trials of
%autocorrelation curves
theta_acmean_std = mean(theta_ac_std);
theta_acsd_std = std(theta_ac_std);
theta_acserr_std = std(theta_ac_std)/sqrt(size(theta_ac_std,1));
gamma_acmean_std = mean(gamma_ac_std);
gamma_acsd_std = std(gamma_ac_std);
gamma_acserr_std = std(gamma_ac_std)/sqrt(size(gamma_ac_std,1));

theta_acmean_sine = mean(theta_ac_sine);
theta_acsd_sine = std(theta_ac_sine);
theta_acserr_sine = std(theta_ac_sine)/sqrt(size(theta_ac_sine,1));
gamma_acmean_sine = mean(gamma_ac_sine);
gamma_acsd_sine = std(gamma_ac_sine);
gamma_acserr_sine = std(gamma_ac_sine)/sqrt(size(gamma_ac_sine,1));

theta_acmean_npt = mean(theta_ac_npt);
theta_acsd_npt = std(theta_ac_npt);
theta_acserr_npt = std(theta_ac_npt)/sqrt(size(theta_ac_npt,1));
gamma_acmean_npt = mean(gamma_ac_npt);
gamma_acsd_npt = std(gamma_ac_npt);
gamma_acserr_npt = std(gamma_ac_npt)/sqrt(size(gamma_ac_npt,1));

theta_acmean_psn = mean(theta_ac_psn);
theta_acsd_psn = std(theta_ac_psn);
theta_acserr_psn = std(theta_ac_psn)/sqrt(size(theta_ac_psn,1));
gamma_acmean_psn = mean(gamma_ac_psn);
gamma_acsd_psn = std(gamma_ac_psn);
gamma_acserr_psn = std(gamma_ac_psn)/sqrt(size(gamma_ac_psn,1));

theta_acmean_rec = mean(theta_ac_rec);
theta_acsd_rec = std(theta_ac_rec);
theta_acserr_rec = std(theta_ac_rec)/sqrt(size(theta_ac_rec,1));
gamma_acmean_rec = mean(gamma_ac_rec);
gamma_acsd_rec = std(gamma_ac_rec);
gamma_acserr_rec = std(gamma_ac_rec)/sqrt(size(gamma_ac_rec,1));


lags = bp_std.autocorr.lags;


cd('/Users/ERCOLE/Documents/Research/Repos/param_analysis/saved_figs');

%% Section for actually making the figures; make sure to add violin plot path first

legend_labels = {'Standard','Sine','Nested Pulse','Poisson','Behavior'};

figure

hold on
% plot(lags,theta_acmean_std)
% plot(lags,theta_acmean_sine)
% plot(lags,theta_acmean_npt)
% plot(lags,theta_acmean_psn)
% plot(lags,theta_acmean_rec)
plot_patch(lags,theta_acmean_sine,theta_acserr_sine,[0 0.4470 0.7410])
plot_patch(lags,theta_acmean_std,theta_acserr_std,[0.8500 0.3250 0.0980])
plot_patch(lags,theta_acmean_psn,theta_acserr_psn,[0.9290 0.6940 0.1250])
plot_patch(lags,theta_acmean_npt,theta_acserr_npt,[0.4940 0.1840 0.5560])
plot_patch(lags,theta_acmean_rec,theta_acserr_rec,'black')
set(gca,'FontSize', 18);
xlabel('Lag (s)')
ylabel('Correlation coefficient')
title('Theta Stim')
legend(legend_labels)
saveas(gcf,'Autocorr_theta.fig')

figure
hold on
% plot(lags,gamma_acmean_std)
% plot(lags,gamma_acmean_sine)
% plot(lags,gamma_acmean_npt)
% plot(lags,gamma_acmean_psn)
% plot(lags,gamma_acmean_rec)
plot_patch(lags,gamma_acmean_sine,gamma_acserr_sine,[0 0.4470 0.7410])
plot_patch(lags,gamma_acmean_std,gamma_acserr_std,[0.8500 0.3250 0.0980])
plot_patch(lags,gamma_acmean_psn,gamma_acserr_psn,[0.9290 0.6940 0.1250])
plot_patch(lags,gamma_acmean_npt,gamma_acserr_npt,[0.4940 0.1840 0.5560])
plot_patch(lags,gamma_acmean_rec,gamma_acserr_rec,'black')
set(gca,'FontSize', 18);
xlabel('Lag (s)')
ylabel('Correlation coefficient')
title('Gamma Stim')
legend(legend_labels)
saveas(gcf,'Autocorr_gamma.fig')

figure
violinplot(theta_phi_vplot,legend_labels)
title('Theta Stim - Phase')
ylabel('Radians')
set(gca,'FontSize', 18);
saveas(gcf,'Phase_vplot_theta.fig')

figure
violinplot(gamma_phi_vplot,legend_labels)
title('Gamma Stim - Phase')
ylabel('Radians')
set(gca,'FontSize', 18);
saveas(gcf,'Phase_vplot_gamma.fig')


function plot_patch(x,mean_tr,dev_tr,col)
    
x_new = [x fliplr(x)];
y_plot = [mean_tr - dev_tr, fliplr(mean_tr + dev_tr)];
patch(x_new, y_plot, col,'FaceAlpha',0.7);

end

