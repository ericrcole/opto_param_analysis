%%
clear all;
%addpath('/Users/ERCOLE/Documents/Research/Data/Behavior')
addpath('/Volumes/OPTO_PROC/OptoStim')

addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/tools/gpml-matlab-v3.6-2015-07-07'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/modeling'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/optogenetic_optimization_project'));

save_figs = true;

% Canonical band ranges
canonical_bands = {'delta','theta','theta2','alpha','beta','lowgamma'};
delta_range = [0.5 4];
theta_range = [4 10];
theta2_range = [4 8];
alpha_range = [8 12];
beta_range = [15 32];
lowgamma_range = [32 50];

stats_grid = {[5,7,11,17,23,29,35,42],
            [5,7,11,17,23,35,42],
            [5,7,11,17,23,29,35,42]};
            
fr_lims = [1200, 1000, 600];
y_star = 0.9*fr_lims;
plot_stars = true;

mean_cb = mean([delta_range; theta_range; theta2_range; alpha_range; beta_range; lowgamma_range],2);


subjects = {'STV003' 'STV008' 'STV009','ARN088'};
states = {'Awake','Autocorr'};
stimuli = {'Standard','Sine','Poisson'};

pvals_all = cell(size(stimuli));
df_all =  cell(size(stimuli));

method = 'pct_change';%'theta_power';%'pct_change';%'theta_power';%'percent_change';%'gamma_power';

%%
% How is the data stored?
% bandpower_data.delta holds the mean delta power over the CA1 and CA3
% channels pre- and post-stim for all trials (N).
% bandpower_data.delta(N,2,2)
% bandpower_data.delta(N,1,1) = mean(mean(delta_pre_CA1))
% bandpower_data.delta(N,2,1) = mean(mean(delta_CA1))
% bandpower_data.delta(N,1,2) = mean(mean(delta_pre_CA3))
% bandpower_data.delta(N,2,2) = mean(mean(delta_CA3))

% bandpower_data.param holds parameters of stimulation
% e.g. bandpower_data.param for 'Standard' holds
% param(k,1) = (amp_v*opto_c(1)+opto_c(2))/(pi*.1^2)  (== Optic Intensity (I))
% param(k,2) = pulse_frequency_pulse
% param(k,3) = stimulation_pulse_width_pulse_a
% e.g. bandpower_data.param for 'NPT' holds
% amp_v = stim_table.stimulation_amplitude_a(k);
% param(k,1) = (amp_v*opto_c(1)+opto_c(2))/(pi*.1^2);
% param(k,2) = stim_table.pulse_frequency_train(k);
% param(k,3) = stim_table.pulse_frequency_pulse(k);
% param(k,4) = stim_table.stimulation_pulse_width_train_a(k);
% param(k,5) = stim_table.stimulation_pulse_width_pulse_a(k);

%plot(bandpower_data.fr_inds(5:end), bandpower_data.stim_psd_CA1{20,1}(5:end))



%% Plotting Loops
f1 = figure('Position',  [100, 100, 1000, 1150]);
% Stimulus
plt = 1;
for i = 1:length(stimuli)
    
    % Awake and Anesthesia states
    for k = 1:length(states)
        stim_freq = [];
        canonical_band = [];
        percent_change = [];
        pct_chg = [];
        stim_prestim = [];
        difference = [];
        pulse_freq_train = [];
        pulse_width_train = [];
        theta_power = [];
        gamma_power = [];
        
        phase_bp = [];
        phase_theta = [];
        phase_gamma = [];
        freq_phase = [];
        
        autocorr_bp = [];
        freq_autocorr = [];
        
        autocorr_rec = cell(0,1);
        
        % Loop through all subjects
        for j = 1:length(subjects)
            load(sprintf('%s_BandpowerPhase_%s_%s.mat',subjects{j},'Awake',stimuli{i}))
            bp_rec = load(sprintf('%s_BandpowerPhase_%s_%s.mat',subjects{j},'Awake','Recording'));
            
            experiment_starts = find(bandpower_data.stim_order == 1);
            % Get baseline band power for all canonical bands as mean, std of prestim
            % bandpower. Average of CA1 and CA3
            experiment_idx = [];
            baseline = [];
            for m = 1:length(experiment_starts)
                experiment_start = experiment_starts(m);
                if m < length(experiment_starts)
                    experiment_end = experiment_starts(m+1);
                    experiment_idx(experiment_start:experiment_end) = m;
                    
                    baseline(1,m) = mean(mean(bandpower_data.delta(experiment_start:experiment_end,1,:),3));
                    baseline(2,m) = mean(mean(bandpower_data.theta(experiment_start:experiment_end,1,:),3));
                    baseline(3,m) = mean(mean(bandpower_data.theta2(experiment_start:experiment_end,1,:),3));
                    baseline(4,m) = mean(mean(bandpower_data.alpha(experiment_start:experiment_end,1,:),3));
                    baseline(5,m) = mean(mean(bandpower_data.beta(experiment_start:experiment_end,1,:),3));
                    baseline(6,m) = mean(mean(bandpower_data.lowgamma(experiment_start:experiment_end,1,:),3));
                else
                    experiment_idx(experiment_start:length(bandpower_data.stim_order)) = m;
                    
                    baseline(1,m) = mean(mean(bandpower_data.delta(experiment_start:end,1,:),3));
                    baseline(2,m) = mean(mean(bandpower_data.theta(experiment_start:end,1,:),3));
                    baseline(3,m) = mean(mean(bandpower_data.theta2(experiment_start:end,1,:),3));
                    baseline(4,m) = mean(mean(bandpower_data.alpha(experiment_start:end,1,:),3));
                    baseline(5,m) = mean(mean(bandpower_data.beta(experiment_start:end,1,:),3));
                    baseline(6,m) = mean(mean(bandpower_data.lowgamma(experiment_start:end,1,:),3));
                end
            end

%             baseline = mean(prestim,2);
%             baseline_std = std(prestim,[],2);
            % Loops through all trials of experiment
            % Find band range a certain trial belongs to and compare the bandpower in
            % that range to baseline
            for trial = 1:length(bandpower_data.stim_times)
                % only do analysis on Intensity > 30
                if bandpower_data.param(trial,1) < 30
                    continue
                end
                exp_idx = experiment_idx(trial);
                if bandpower_data.param(trial,2) > delta_range(1) && ...
                        bandpower_data.param(trial,2) < delta_range(2)
                    stim_freq = [stim_freq bandpower_data.param(trial,2)];
                    if strcmp(stimuli{i}, 'NPT')
                        pulse_freq_train = [pulse_freq_train bandpower_data.param(trial,3)];
                        pulse_width_train = [pulse_width_train bandpower_data.param(trial,4)];
                    end
                    theta = mean(bandpower_data.theta(trial,2,:),3);
                    theta_prestim = mean(bandpower_data.theta(trial,1,:),3);
                    gamma = mean(bandpower_data.lowgamma(trial,2,:),3);
                    gamma_prestim = mean(bandpower_data.lowgamma(trial,1,:),3);
                    %theta_power = [theta_power (theta/baseline(2,exp_idx)-1)*100];
                    %gamma_power = [gamma_power (gamma/baseline(6,exp_idx)-1)*100];
                    canonical_band = [canonical_band mean(delta_range)];
                    bandpower = mean(bandpower_data.delta(trial,2,:),3);
                    percent_change = [percent_change (bandpower/baseline(1,exp_idx) - 1)*100];
                    bandpower_prestim = mean(bandpower_data.delta(trial,1,:),3);
                    pct_chg = [pct_chg (bandpower/bandpower_prestim - 1)*100];
                    stim_prestim = [stim_prestim (bandpower-bandpower_prestim)/baseline(1)*100];
                    difference = [difference bandpower-baseline(1)];
                    theta_power = [theta_power (theta/theta_prestim - 1)*100];
                    gamma_power = [gamma_power (gamma/gamma_prestim - 1)*100];
                    
                    
                    
                elseif bandpower_data.param(trial,2) > theta_range(1) && ...
                        bandpower_data.param(trial,2) < theta_range(2)
                    stim_freq = [stim_freq bandpower_data.param(trial,2)];
                    if strcmp(stimuli{i}, 'NPT')
                        pulse_freq_train = [pulse_freq_train bandpower_data.param(trial,3)];
                        pulse_width_train = [pulse_width_train bandpower_data.param(trial,4)];
                    end
                    theta = mean(bandpower_data.theta(trial,2,:),3);
                    theta_prestim = mean(bandpower_data.theta(trial,1,:),3);
                    gamma = mean(bandpower_data.lowgamma(trial,2,:),3);
                    gamma_prestim = mean(bandpower_data.lowgamma(trial,1,:),3);
%                     theta_power = [theta_power (mean(bandpower_data.theta(trial,2,:),3)/baseline(2,exp_idx)-1)*100];
%                     gamma_power = [gamma_power (mean(bandpower_data.lowgamma(trial,2,:),3)/baseline(6,exp_idx)-1)*100];
                    canonical_band = [canonical_band mean(theta_range)];
                    bandpower = mean(bandpower_data.theta(trial,2,:),3);
                    percent_change = [percent_change (bandpower/baseline(2,exp_idx) - 1)*100];
                    bandpower_prestim = mean(bandpower_data.theta(trial,1,:),3);
                    pct_chg = [pct_chg (bandpower/bandpower_prestim - 1)*100];
                    stim_prestim = [stim_prestim (bandpower-bandpower_prestim)/baseline(2)*100];
                    difference = [difference bandpower-baseline(2)];
                    theta_power = [theta_power (theta/theta_prestim - 1)*100];
                    gamma_power = [gamma_power (gamma/gamma_prestim - 1)*100];
                elseif bandpower_data.param(trial,2) > theta2_range(1) && ...
                        bandpower_data.param(trial,2) < theta2_range(2)
                    stim_freq = [stim_freq bandpower_data.param(trial,2)];
                    if strcmp(stimuli{i}, 'NPT')
                        pulse_freq_train = [pulse_freq_train bandpower_data.param(trial,3)];
                        pulse_width_train = [pulse_width_train bandpower_data.param(trial,4)];
                    end
                    theta = mean(bandpower_data.theta(trial,2,:),3);
                    theta_prestim = mean(bandpower_data.theta(trial,1,:),3);
                    gamma = mean(bandpower_data.lowgamma(trial,2,:),3);
                    gamma_prestim = mean(bandpower_data.lowgamma(trial,1,:),3);
%                     theta_power = [theta_power (mean(bandpower_data.theta(trial,2,:),3)/baseline(2,exp_idx)-1)*100];
%                     gamma_power = [gamma_power (mean(bandpower_data.lowgamma(trial,2,:),3)/baseline(6,exp_idx)-1)*100];
                    canonical_band = [canonical_band mean(theta2_range)];
                    bandpower = mean(bandpower_data.theta2(trial,2,:),3);
                    percent_change = [percent_change (bandpower/baseline(3,exp_idx) - 1)*100];
                    bandpower_prestim = mean(bandpower_data.theta2(trial,1,:),3);
                    pct_chg = [pct_chg (bandpower/bandpower_prestim - 1)*100];
                    stim_prestim = [stim_prestim (bandpower-bandpower_prestim)/baseline(3)*100];
                    difference = [difference bandpower-baseline(3)];
                    theta_power = [theta_power (theta/theta_prestim - 1)*100];
                    gamma_power = [gamma_power (gamma/gamma_prestim - 1)*100];
                elseif bandpower_data.param(trial,2) > alpha_range(1) && ...
                        bandpower_data.param(trial,2) < alpha_range(2)
                    stim_freq = [stim_freq bandpower_data.param(trial,2)];
                    if strcmp(stimuli{i}, 'NPT')
                        pulse_freq_train = [pulse_freq_train bandpower_data.param(trial,3)];
                        pulse_width_train = [pulse_width_train bandpower_data.param(trial,4)];
                    end
                    theta = mean(bandpower_data.theta(trial,2,:),3);
                    theta_prestim = mean(bandpower_data.theta(trial,1,:),3);
                    gamma = mean(bandpower_data.lowgamma(trial,2,:),3);
                    gamma_prestim = mean(bandpower_data.lowgamma(trial,1,:),3);
%                     theta_power = [theta_power (mean(bandpower_data.theta(trial,2,:),3)/baseline(2,exp_idx)-1)*100];
%                     gamma_power = [gamma_power (mean(bandpower_data.lowgamma(trial,2,:),3)/baseline(6,exp_idx)-1)*100];
                    canonical_band = [canonical_band mean(alpha_range)];
                    bandpower = mean(bandpower_data.alpha(trial,2,:),3);
                    percent_change = [percent_change (bandpower/baseline(4,exp_idx) - 1)*100];
                    bandpower_prestim = mean(bandpower_data.alpha(trial,1,:),3);
                    pct_chg = [pct_chg (bandpower/bandpower_prestim - 1)*100];
                    stim_prestim = [stim_prestim (bandpower-bandpower_prestim)/baseline(4)*100];
                    difference = [difference bandpower-baseline(4)];
                    theta_power = [theta_power (theta/theta_prestim - 1)*100];
                    gamma_power = [gamma_power (gamma/gamma_prestim - 1)*100];
                elseif bandpower_data.param(trial,2) > beta_range(1) && ...
                        bandpower_data.param(trial,2) < beta_range(2)
                    stim_freq = [stim_freq bandpower_data.param(trial,2)];
                    if strcmp(stimuli{i}, 'NPT')
                        pulse_freq_train = [pulse_freq_train bandpower_data.param(trial,3)];
                        pulse_width_train = [pulse_width_train bandpower_data.param(trial,4)];
                    end
                    theta = mean(bandpower_data.theta(trial,2,:),3);
                    theta_prestim = mean(bandpower_data.theta(trial,1,:),3);
                    gamma = mean(bandpower_data.lowgamma(trial,2,:),3);
                    gamma_prestim = mean(bandpower_data.lowgamma(trial,1,:),3);
%                     theta_power = [theta_power (mean(bandpower_data.theta(trial,2,:),3)/baseline(2,exp_idx)-1)*100];
%                     gamma_power = [gamma_power (mean(bandpower_data.lowgamma(trial,2,:),3)/baseline(6,exp_idx)-1)*100];
                    canonical_band = [canonical_band mean(beta_range)];
                    bandpower = mean(bandpower_data.beta(trial,2,:),3);
                    percent_change = [percent_change (bandpower/baseline(5,exp_idx) - 1)*100];
                    bandpower_prestim = mean(bandpower_data.beta(trial,1,:),3);
                    pct_chg = [pct_chg (bandpower/bandpower_prestim - 1)*100];
                    stim_prestim = [stim_prestim (bandpower-bandpower_prestim)/baseline(5)*100];
                    difference = [difference bandpower-baseline(5)];
                    theta_power = [theta_power (theta/theta_prestim - 1)*100];
                    gamma_power = [gamma_power (gamma/gamma_prestim - 1)*100];
                elseif bandpower_data.param(trial,2) > lowgamma_range(1) && ...
                        bandpower_data.param(trial,2) < lowgamma_range(2)
                    stim_freq = [stim_freq bandpower_data.param(trial,2)];
                    if strcmp(stimuli{i}, 'NPT')
                        pulse_freq_train = [pulse_freq_train bandpower_data.param(trial,3)];
                        pulse_width_train = [pulse_width_train bandpower_data.param(trial,4)];
                    end
                    theta = mean(bandpower_data.theta(trial,2,:),3);
                    theta_prestim = mean(bandpower_data.theta(trial,1,:),3);
                    gamma = mean(bandpower_data.lowgamma(trial,2,:),3);
                    gamma_prestim = mean(bandpower_data.lowgamma(trial,1,:),3);
%                     theta_power = [theta_power (mean(bandpower_data.theta(trial,2,:),3)/baseline(2,exp_idx)-1)*100];
%                     gamma_power = [gamma_power (mean(bandpower_data.lowgamma(trial,2,:),3)/baseline(6,exp_idx)-1)*100];
                    canonical_band = [canonical_band mean(lowgamma_range)];
                    bandpower = mean(bandpower_data.lowgamma(trial,2,:),3);
                    percent_change = [percent_change (bandpower/baseline(6,exp_idx) - 1)*100];
                    bandpower_prestim = mean(bandpower_data.lowgamma(trial,1,:),3);
                    pct_chg = [pct_chg (bandpower/bandpower_prestim - 1)*100];
                    stim_prestim = [stim_prestim (bandpower-bandpower_prestim)/baseline(6)*100];
                    difference = [difference bandpower-baseline(6)];
                    theta_power = [theta_power (theta/theta_prestim - 1)*100];
                    gamma_power = [gamma_power (gamma/gamma_prestim - 1)*100];
                end
            end
            %stim_freq = stim frequencies
            amp_inds = bandpower_data.param(:,1)>=30;
            [phase_temp,freqs_temp] = get_stim_phase(bandpower_data.param(amp_inds,2),bandpower_data.phi.stim(amp_inds),bandpower_data.phi.fr_inds);
            phase_bp = [phase_bp; phase_temp];
            freq_phase = [freq_phase; freqs_temp];
            
            [autocorr_temp,freqs_temp] = get_stim_autocorr(bandpower_data.param(amp_inds,2),bandpower_data.autocorr.stim(amp_inds),bandpower_data.autocorr.lags);
            autocorr_bp = [autocorr_bp; autocorr_temp];
            freq_autocorr = [freq_autocorr; freqs_temp];
            
            autocorr_rec = [autocorr_rec; bp_rec.bandpower_data.autocorr.stim];
            
        end
        
        stim_freq = stim_freq';
        canonical_band = canonical_band';
        percent_change = percent_change';
        pct_chg = pct_chg';
        stim_prestim = stim_prestim';
        pulse_freq_train = pulse_freq_train';
        pulse_freq_pulse = stim_freq;
        theta_power = theta_power';
        gamma_power = gamma_power';
        
        autocorr_rec_mean = mean(cell2mat(cellfun(@transpose,autocorr_rec,'UniformOutput',false)));
        autocorr_rec_sd = std(cell2mat(cellfun(@transpose,autocorr_rec,'UniformOutput',false)));
        autocorr_rec_serr = autocorr_rec_sd./sqrt(size(autocorr_rec,1));

        % Fit Gaussian Processes
        if strcmp(stimuli{i}, 'NPT')
            predictors = [pulse_freq_train pulse_freq_pulse];
            idxs_high = pulse_width_train > 0.05;
            idxs_low = ~(pulse_width_train > 0.05);
            
            predictors_high = predictors(idxs_high,:);
            predictors_low = predictors(idxs_low,:);
            if strcmp(method,'percent_change')
                percent_change_high = percent_change(idxs_high);
                percent_change_low = percent_change(idxs_low);
                scatterpts_high = percent_change(idxs_high);
                scatterpts_low = percent_change(idxs_low);
            elseif strcmp(method,'pct_change')
                percent_change_high = pct_chg(idxs_high);
                percent_change_low = pct_chg(idxs_low);
                scatterpts_high = pct_chg(idxs_high);
                scatterpts_low = pct_chg(idxs_low);
            elseif strcmp(method,'gamma_power')
                percent_change_high = gamma_power(idxs_high);
                percent_change_low = gamma_power(idxs_low);
                scatterpts_high = gamma_power(idxs_high);
                scatterpts_low = gamma_power(idxs_low);
            elseif strcmp(method,'theta_power')
                percent_change_high = theta_power(idxs_high);
                percent_change_low = theta_power(idxs_low);
                scatterpts_high = theta_power(idxs_high);
                scatterpts_low = theta_power(idxs_low);
%                 idxs = find(isnan(percent_change_high));
%                 idxs1 = find(isinf(percent_change_low));
%                 predictors([idxs; idxs1]) = [];
%                 pct_chg([idxs; idxs1]) = [];
%                 gp1.initialize_data(predictors,pct_chg,min(predictors),max(predictors));
%                 scatterpts = pct_chg;
            end
            gp1 = gp_object();
            gp2 = gp_object();
            
            gp1.initialize_data(predictors_high,percent_change_high,min(predictors_high),max(predictors_high));
            gp2.initialize_data(predictors_low,percent_change_low,min(predictors_low),max(predictors_low));
            
            for it = 1:20
                % Remove outlier scatterpoints
                idxs = find(abs((scatterpts_high - mean(scatterpts_high))/std(scatterpts_high)) > 3.5);
                scatterpts_high(idxs) = [];
                %stim_freq(idxs) = [];
                predictors_high(idxs,:) = [];
                gp1.initialize_data(predictors_high,scatterpts_high,min(predictors_high),max(predictors_high));

                idxs = find(abs((scatterpts_low - mean(scatterpts_low))/std(scatterpts_low)) > 3.5);
                scatterpts_low(idxs) = [];
                %stim_freq(idxs) = [];
                predictors_low(idxs,:) = [];
                gp2.initialize_data(predictors_low,scatterpts_low,min(predictors_low),max(predictors_low));
            end
        elseif strcmp(states{k},'Phase')
            phase_bp(phase_bp<0) = phase_bp(phase_bp<0)+pi;
            predictors = freq_phase;
            pct_change = phase_bp;
            percent_change = phase_bp;
            scatterpts = phase_bp;
            stim_freq = freq_phase;
            
            gp1.initialize_data(predictors,pct_chg,min(predictors),max(predictors));
            for it = 1:20
                % Remove outlier scatterpoints
                idxs = find(abs((scatterpts - mean(scatterpts))/std(scatterpts)) > 3.5);
                scatterpts(idxs) = [];
                stim_freq(idxs) = [];
                predictors(idxs,:) = [];
                gp1.initialize_data(predictors,scatterpts,min(predictors),max(predictors));
            end
        elseif strcmp(states{k},'Autocorr')
            
            predictors = freq_phase;
            pct_change = autocorr_bp;
            percent_change = autocorr_bp;
            scatterpts = autocorr_bp;
            stim_freq = freq_autocorr;
            
            gp1.initialize_data(predictors,pct_chg,min(predictors),max(predictors));
            for it = 1:20
                % Remove outlier scatterpoints
                idxs = find(abs((scatterpts - mean(scatterpts))/std(scatterpts)) > 3.5);
                scatterpts(idxs) = [];
                stim_freq(idxs) = [];
                predictors(idxs,:) = [];
                gp1.initialize_data(predictors,scatterpts,min(predictors),max(predictors));
            end
        else
            predictors = stim_freq;
            gp1 = gp_object();
            if strcmp(method,'percent_change')
                gp1.initialize_data(predictors,percent_change,min(predictors),max(predictors));
                scatterpts = percent_change;
            elseif strcmp(method,'pct_change')
                idxs = find(isnan(pct_chg));
                idxs1 = find(isinf(pct_chg));
                stim_freq([idxs; idxs1]) = [];
                predictors([idxs; idxs1]) = [];
                pct_chg([idxs; idxs1]) = [];
                gp1.initialize_data(predictors,pct_chg,min(predictors),max(predictors));
                scatterpts = pct_chg;
            elseif strcmp(method,'stim_prestim')
                gp1.initialize_data(predictors,stim_prestim,min(predictors),max(predictors));
                scatterpts = stim_prestim;
            elseif strcmp(method,'difference')
                gp1.initialize_data(predictors,difference,min(predictors),max(predictors));
                scatterpts = difference;
            elseif strcmp(method,'gamma_power')
                idxs = find(isnan(gamma_power));
                idxs1 = find(isinf(gamma_power));
                stim_freq([idxs; idxs1]) = [];
                predictors([idxs; idxs1]) = [];
                gamma_power([idxs; idxs1]) = [];
                gp1.initialize_data(predictors,gamma_power,min(predictors),max(predictors));
                scatterpts = gamma_power;
            elseif strcmp(method,'theta_power')
                idxs = find(isnan(theta_power));
                idxs1 = find(isinf(theta_power));
                stim_freq([idxs; idxs1]) = [];
                predictors([idxs; idxs1]) = [];
                theta_power([idxs; idxs1]) = [];
                gp1.initialize_data(predictors,theta_power,min(predictors),max(predictors));
                scatterpts = theta_power;
            end
            for it = 1:20
                % Remove outlier scatterpoints
                idxs = find(abs((scatterpts - mean(scatterpts))/std(scatterpts)) > 3.5);
                scatterpts(idxs) = [];
                stim_freq(idxs) = [];
                predictors(idxs,:) = [];
                gp1.initialize_data(predictors,scatterpts,min(predictors),max(predictors));
            end
%             % Remove outlier scatterpoints
%             idxs = find(abs((scatterpts - mean(scatterpts))/std(scatterpts)) > 2);
%             scatterpts(idxs) = [];
%             stim_freq(idxs) = [];
%             predictors(idxs,:) = [];
%             gp1.initialize_data(predictors,scatterpts,min(predictors),max(predictors));
        end
        
        if any(contains(stimuli,'NPT'))
            subplot(length(stimuli)+1,length(states),plt)
        else
            subplot(length(stimuli),length(states),plt)
        end
        
        % STATS BLOCK
        if strcmp(states{k},'Awake')
            pvals = nan(size(stats_grid{i}));
            df_temp = nan(size(stats_grid{i}));
            for pind = 1:length(stats_grid{i})
                gridinds = find(predictors == stats_grid{i}(pind));
                gridvals = scatterpts(gridinds);
                [pvals(pind),h,stats] = signrank(gridvals);
                df_temp(pind) = length(gridinds);
            end
            pvals_all{i} = pvals;
            df_all{i} = df_temp;
        end
        
        if strcmp(stimuli{i}, 'NPT')
            %figure()
            plot_surface(gp1)
            ylabel('Pulse Train Frequency (Hz)')
            xlabel('Pulse Frequency (Hz)')
            title([states{k} ', ' stimuli{i} ' Train Width > 0.05'],'FontSize',14)
            if strcmp(method,'gamma_power')
                zlabel('Change in Gamma (%)','FontSize',14)
            elseif strcmp(method,'theta_power')
                zlabel('Change in Theta (%)','FontSize',14)
            else
                zlabel('Change from baseline (%)','FontSize',14)
            end
            if any(contains(stimuli,'NPT'))
                subplot(length(stimuli)+1,length(states),plt)
            else
                subplot(length(stimuli),length(states),plt+2)
            end
            
            %figure()
            plot_surface(gp2)
            ylabel('Pulse Train Frequency (Hz)')
            xlabel('Pulse Frequency (Hz)')
            if strcmp(method,'gamma_power')
                zlabel('Change in Gamma (%)','FontSize',14)
            elseif strcmp(method,'theta_power')
                zlabel('Change in Theta (%)','FontSize',14)
            else
                zlabel('Change from baseline (%)','FontSize',14)
            end
            title([states{k} ', ' stimuli{i} ' Train Width <= 0.05'],'FontSize',14)
        else
            %PLOT SCATTER AND GP FIT
            [y_pred, y_std] = plot_fit(gp1);
            hold on
            scatter(stim_freq,scatterpts,5,'k','filled')
            if plot_stars
                for kk = 1:length(pvals)
                   if pvals(kk) < 0.0001
                       text(stats_grid{i}(kk),y_star(i),'***','FontSize',28,'HorizontalAlignment','center')
                   elseif pvals(kk) < 0.001
                       text(stats_grid{i}(kk),y_star(i),'**','FontSize',28,'HorizontalAlignment','center')
                   elseif pvals(kk) < 0.01
                       text(stats_grid{i}(kk),y_star(i),'*','FontSize',28,'HorizontalAlignment','center')
                   else
                       %text(stats_grid{i}(kk),y_star(i),'n','FontSize',14,'HorizontalAlignment','center')
                   end
                end
            end
            hold off
            if i == 3
                xlabel('Stimulation Frequency (Hz)','FontSize',22)
            end
            xlim([-.5 50.5])
            
            if strcmp(states{k},'Phase')
                ylim([-.1 3.2])
                if strcmp(method,'gamma_power')
                    ylabel('Change in Gamma (%)','FontSize',14)
                elseif strcmp(method,'theta_power')
                    ylabel('Change in Theta (%)','FontSize',14)
                else
                    ylabel('Phase (radians)','FontSize',14)
                end
            elseif strcmp(states{k},'Autocorr')
                ylim([-.05 1.1])
                hold on
                lag_ind1 = find(bandpower_data.autocorr.lags >= 1/50,1);
                lag_ind2 = find(bandpower_data.autocorr.lags >= 1/2,1);
                lags = 1./(bandpower_data.autocorr.lags(lag_ind1:lag_ind2));
                
                %%%plot mean of autocorr
                %p1 = plot(lags,autocorr_rec_mean(lag_ind1:lag_ind2),'red','LineWidth',1.5);
                plot_patch(lags,autocorr_rec_mean(lag_ind1:lag_ind2),autocorr_rec_serr(lag_ind1:lag_ind2),'red');
                if strcmp(method,'gamma_power')
                    ylabel('Change in Gamma (%)','FontSize',18)
                elseif strcmp(method,'theta_power')
                    ylabel('Change in Theta (%)','FontSize',18)
                else
                    ylabel('Correlation coefficient','FontSize',18)
                end
                if i == 1
                    title('Autocorrelation','FontSize',20)
                end
            elseif strcmp(states{k},'Awake')
                ylim([-200 fr_lims(i)])
                if strcmp(method,'gamma_power')
                    ylabel('Change in Gamma (%)','FontSize',18)
                elseif strcmp(method,'theta_power')
                    ylabel('Change in Theta (%)','FontSize',18)
                else
                    ylabel({stimuli{i},' ','Percent change (%)'},'FontSize',22)
                end
                if i == 1
                    title('Bandpower entrainment','FontSize',28)
                end
            end
            
            %title([states{k} ', ' stimuli{i}],'FontSize',14)
        end
        set(gca,'fontsize',15)
        % Remove outliers for plotting only
        %ylim([-100 1000])
        
        
        plt = plt+1;

    end
end

annotation('textbox',[.03 .88 .1 .1],'String','a','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 28);
annotation('textbox',[.03 .58 .1 .1],'String','b','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 28);
annotation('textbox',[.03 .285 .1 .1],'String','c','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 28);

if save_figs
    print(gcf,'fig3.png','-dpng','-r300')
%saveas(gcf,'fig3.png')
%saveas(gcf,'fig3.fig')
end

function [phase_out,new_freqs] = get_stim_phase(stim_freq, bp_phase, fr_inds)
    
    val_inds = find(~cellfun(@isempty,bp_phase));
    stim_freq = stim_freq(val_inds);
    bp_phase = bp_phase(val_inds);
    phase_out = zeros(size(stim_freq));
    
    exc_inds = zeros(size(phase_out));
    
    for k = 1:length(phase_out)
        stim_ind = find(fr_inds>=stim_freq(k),1);
        if stim_ind>length(bp_phase{k})
            exc_inds(k) = 1;
            continue
        end
        if stim_freq(k) == fr_inds(stim_ind)
           phase_out(k) = bp_phase{k}(stim_ind);
        else
           phase_out(k) = (bp_phase{k}(stim_ind)+bp_phase{k}(stim_ind-1))/2;
        end
    end
    phase_out(logical(exc_inds)) = [];
    stim_freq(logical(exc_inds)) = [];
    new_freqs = stim_freq;

end


function [phase_out,new_freqs] = get_stim_autocorr(stim_freq, bp_phase, fr_inds)
    
    val_inds = find(~cellfun(@isempty,bp_phase));
    stim_freq = stim_freq(val_inds);
    bp_phase = bp_phase(val_inds);
    phase_out = zeros(size(stim_freq));
    
    exc_inds = zeros(size(phase_out));
    
    for k = 1:length(phase_out)
        stim_ind = find(fr_inds>=(1/stim_freq(k)),1);
        if stim_ind>length(bp_phase{k})
            exc_inds(k) = 1;
            continue
        end
        if stim_freq(k) == fr_inds(stim_ind)
           phase_out(k) = bp_phase{k}(stim_ind);
        else
           phase_out(k) = (bp_phase{k}(stim_ind)+bp_phase{k}(stim_ind-1))/2;
        end
    end
    phase_out(logical(exc_inds)) = [];
    stim_freq(logical(exc_inds)) = [];
    new_freqs = stim_freq;

end

function plot_patch(x,mean_tr,dev_tr,col)
    
x_new = [x fliplr(x)];
y_plot = [mean_tr - dev_tr, fliplr(mean_tr + dev_tr)];
patch(x_new, y_plot, col,'FaceAlpha',0.4);

end