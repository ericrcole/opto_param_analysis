addpath('/Users/ERCOLE/Documents/Research/Data/Behavior/Phase')

addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/tools/gpml-matlab-v3.6-2015-07-07'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/modeling'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/optogenetic_optimization_project'));
addpath('Violinplot-Matlab')

subjects = {'ARN088' 'STV003' 'STV008' 'STV009'};
states = {'Awake','Anesthesia'};
stimuli = {'Standard','Sine','NPT'};

kk = 2;

load(sprintf('%s_BandpowerPhase_Awake_Standard.mat',subjects{kk}));
bp_std = bandpower_data;
load(sprintf('%s_BandpowerPhase_Awake_Sine.mat',subjects{kk}));
bp_sine = bandpower_data;
load(sprintf('%s_BandpowerPhase_Awake_NPT.mat',subjects{kk}));
bp_npt = bandpower_data;
load(sprintf('%s_BandpowerPhase_Awake_Poisson.mat',subjects{kk}));
bp_psn = bandpower_data;
load(sprintf('%s_BandpowerPhase_Awake_Sine_2N.mat',subjects{kk}));
bp_sine2n = bandpower_data;

labels = {'Pre-stim','Stim'};

theta_inds = find((bp_std.param(:,2)==7));
gamma_inds = find((bp_std.param(:,2)==42));

figure()
%histograms of pre-stim vs during stim phase calc
subplot(2,2,1)
hold on
histogram(mean(bp_std.phi.theta(theta_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_std.phi.theta(theta_inds,2,:),3),20,'Normalization','pdf');
title('Theta'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
subplot(2,2,2)
hold on
histogram(mean(bp_std.phi.lowgamma(gamma_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_std.phi.lowgamma(gamma_inds,2,:),3),20,'Normalization','pdf');
title('Low Gamma'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
suptitle('Phase Entrainment: Standard')

subplot(2,2,3)

violinplot([mean(bp_std.phi.theta(theta_inds,1,:),3),mean(bp_std.phi.theta(theta_inds,2,:),3)],labels);
ylabel('Radians')

subplot(2,2,4)
violinplot([mean(bp_std.phi.lowgamma(gamma_inds,1,:),3),mean(bp_std.phi.lowgamma(gamma_inds,2,:),3)],labels);
ylabel('Radians')

theta_inds = find((bp_sine.param(:,2)==7));
gamma_inds = find((bp_sine.param(:,2)==42));

figure()
subplot(2,2,1)
hold on
histogram(mean(bp_sine.phi.theta(theta_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_sine.phi.theta(theta_inds,2,:),3),20,'Normalization','pdf');
title('Theta'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
subplot(2,2,2)
hold on
histogram(mean(bp_sine.phi.lowgamma(gamma_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_sine.phi.lowgamma(gamma_inds,2,:),3),20,'Normalization','pdf');
title('Low Gamma'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
suptitle('Phase Entrainment: Sine')

subplot(2,2,3)
violinplot([mean(bp_sine.phi.theta(theta_inds,1,:),3),mean(bp_sine.phi.theta(theta_inds,2,:),3)],labels);
ylabel('Radians')

subplot(2,2,4)
violinplot([mean(bp_sine.phi.lowgamma(gamma_inds,1,:),3),mean(bp_sine.phi.lowgamma(gamma_inds,2,:),3)],labels);
ylabel('Radians')

theta_inds = find((bp_psn.param(:,2)==7));
gamma_inds = find((bp_psn.param(:,2)==42));

figure()
subplot(2,2,1)
hold on
histogram(mean(bp_psn.phi.theta(theta_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_psn.phi.theta(theta_inds,2,:),3),20,'Normalization','pdf');
title('Theta'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
subplot(2,2,2)
hold on
histogram(mean(bp_psn.phi.lowgamma(gamma_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_psn.phi.lowgamma(gamma_inds,2,:),3),20,'Normalization','pdf');
title('Low Gamma'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
suptitle('Phase Entrainment: Poisson')

subplot(2,2,3)
violinplot([mean(bp_psn.phi.theta(theta_inds,1,:),3),mean(bp_psn.phi.theta(theta_inds,2,:),3)],labels);
ylabel('Radians')

subplot(2,2,4)
violinplot([mean(bp_psn.phi.lowgamma(gamma_inds,1,:),3),mean(bp_psn.phi.lowgamma(gamma_inds,2,:),3)],labels);
ylabel('Radians')

theta_inds = find(and(bp_npt.param(:,2)<9,bp_npt.param(:,2)>5));
gamma_inds = find(and(bp_npt.param(:,3)<50,bp_npt.param(:,3)>30));

figure()
subplot(2,2,1)
hold on
histogram(mean(bp_npt.phi.theta(theta_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_npt.phi.theta(theta_inds,2,:),3),20,'Normalization','pdf');
title('Theta'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
subplot(2,2,2)
hold on
histogram(mean(bp_npt.phi.lowgamma(gamma_inds,1,:),3),20,'Normalization','pdf');
histogram(mean(bp_npt.phi.lowgamma(gamma_inds,2,:),3),20,'Normalization','pdf');
title('Low Gamma'); ylabel('Probability Density');
legend({'Pre-Stim','Stim'})
suptitle('Phase Entrainment: Nested Pulse')

subplot(2,2,3)
violinplot([mean(bp_npt.phi.theta(theta_inds,1,:),3),mean(bp_npt.phi.theta(theta_inds,2,:),3)],labels);
ylabel('Radians')

subplot(2,2,4)
violinplot([mean(bp_npt.phi.lowgamma(gamma_inds,1,:),3),mean(bp_npt.phi.lowgamma(gamma_inds,2,:),3)],labels);
ylabel('Radians')

