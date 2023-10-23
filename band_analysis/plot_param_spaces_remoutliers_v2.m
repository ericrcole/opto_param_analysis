%plot_stim_standard.m
%
% Plot theta and gamma vs. stimulation parameters.
% Fit GP model and visualize surface plot of response
%
% V2: added visualization of phase for theta/gamma

addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/tools/gpml-matlab-v3.6-2015-07-07'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/modeling'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/optogenetic_optimization_project'));

save = 0;

animal = {'ARN088','STV008','STV003','STV009'};

for kk = 1:length(animal)
cd('/Users/ERCOLE/OneDrive - Emory University/GrossLab_OptoStimData/Processed/Behavior');

load(sprintf('%s_Bandpower_Awake_Standard.mat',animal{kk}));
bp_standard = bandpower_data;
load(sprintf('%s_Bandpower_Awake_Sine.mat',animal{kk}));
bp_sine = bandpower_data;
load(sprintf('%s_Bandpower_Awake_Sine_2N.mat',animal{kk}));
bp_sine2n = bandpower_data;
load(sprintf('%s_Bandpower_Awake_NPT.mat',animal{kk}));
bp_npt = bandpower_data;
load(sprintf('%s_Bandpower_Awake_Poisson.mat',animal{kk}));
bp_poisson = bandpower_data;

cd('/Users/ERCOLE/OneDrive - Emory University/GrossLab_OptoStimData/Code/Repos/param_analysis')

gamma_norm_std = normalize_bp(bp_standard,'gamma');
theta_norm_std = normalize_bp(bp_standard,'theta');
gamma_norm_sine = normalize_bp(bp_sine,'gamma');
theta_norm_sine = normalize_bp(bp_sine,'theta');
gamma_norm_npt = normalize_bp(bp_npt,'gamma');
theta_norm_npt = normalize_bp(bp_npt,'theta');
gamma_norm_psn = normalize_bp(bp_poisson,'gamma');
theta_norm_psn = normalize_bp(bp_poisson,'theta');
gamma_norm_sine2n = normalize_bp(bp_sine2n,'gamma');
theta_norm_sine2n = normalize_bp(bp_sine2n,'theta');

outl_temp = isoutlier(gamma_norm_std,'median','ThresholdFactor',4);
gp_gamma_std = gp_object();
gp_gamma_std.initialize_data(bp_standard.param(~outl_temp,:),gamma_norm_std(~outl_temp), min(bp_standard.param), max(bp_standard.param));
gamma_norm_std(outl_temp) = [];
outl_temp = isoutlier(theta_norm_std,'median','ThresholdFactor',4);
gp_theta_std = gp_object();
gp_theta_std.initialize_data(bp_standard.param(~outl_temp,:),theta_norm_std(~outl_temp), min(bp_standard.param), max(bp_standard.param));
theta_norm_std(outl_temp) = [];

outl_temp = isoutlier(gamma_norm_sine,'median','ThresholdFactor',4);
gp_gamma_sine = gp_object();
gp_gamma_sine.initialize_data(bp_sine.param(~outl_temp,:),gamma_norm_sine(~outl_temp), min(bp_sine.param), max(bp_sine.param));
gamma_norm_sine(outl_temp) = [];
outl_temp = isoutlier(theta_norm_sine,'median','ThresholdFactor',4);
gp_theta_sine = gp_object();
gp_theta_sine.initialize_data(bp_sine.param(~outl_temp,:),theta_norm_sine(~outl_temp), min(bp_sine.param), max(bp_sine.param));
theta_norm_sine(outl_temp) = [];

outl_temp = isoutlier(gamma_norm_npt,'median','ThresholdFactor',4);
gp_gamma_npt = gp_object();
gp_gamma_npt.initialize_data(bp_npt.param(~outl_temp,:),gamma_norm_npt(~outl_temp), min(bp_npt.param), max(bp_npt.param));
gamma_norm_npt(outl_temp) = [];
outl_temp = isoutlier(theta_norm_npt,'median','ThresholdFactor',4);
gp_theta_npt = gp_object();
gp_theta_npt.initialize_data(bp_npt.param(~outl_temp,:),theta_norm_npt(~outl_temp), min(bp_npt.param), max(bp_npt.param));
theta_norm_npt(outl_temp) = [];

outl_temp = isoutlier(gamma_norm_psn,'median','ThresholdFactor',4);
gp_gamma_psn = gp_object();
gp_gamma_psn.initialize_data(bp_poisson.param(~outl_temp,:),gamma_norm_psn(~outl_temp), min(bp_poisson.param), max(bp_poisson.param));
gamma_norm_psn(outl_temp) = [];
outl_temp = isoutlier(theta_norm_psn,'median','ThresholdFactor',4);
gp_theta_psn = gp_object();
gp_theta_psn.initialize_data(bp_poisson.param(~outl_temp,:),theta_norm_psn(~outl_temp), min(bp_poisson.param), max(bp_poisson.param));
theta_norm_psn(outl_temp) = [];

outl_temp = isoutlier(gamma_norm_sine2n,'median','ThresholdFactor',4);
gp_gamma_sine2n = gp_object();
gp_gamma_sine2n.initialize_data(bp_sine2n.param(~outl_temp,:),gamma_norm_sine2n(~outl_temp), min(bp_sine2n.param), max(bp_sine2n.param));
gamma_norm_sine2n(outl_temp) = [];
outl_temp = isoutlier(theta_norm_sine2n,'median','ThresholdFactor',4);
gp_theta_sine2n = gp_object();
gp_theta_sine2n.initialize_data(bp_sine2n.param(~outl_temp,:),theta_norm_sine2n(~outl_temp), min(bp_sine2n.param), max(bp_sine2n.param));
theta_norm_sine2n(outl_temp) = [];

f1 = figure;
subplot(1,2,1)
hold on
%plot3(bp_standard.param(:,1),bp_standard.param(:,2),gamma_norm_std,'bo');
plot_surface(gp_gamma_std);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Gamma Power')
%view([-10,-10,8])
view(45,15)
subplot(1,2,2)
hold on
%plot3(bp_standard.param(:,1),bp_standard.param(:,2),theta_norm_std,'bo');
plot_surface(gp_theta_std);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Theta Power')
%view([-10,-10,8])
view(45,15)
f1.Position = [100 100 1000 600];
suptitle('Standard Pulse')
if save == 1
    saveas(gcf,sprintf('./saved_figs/%s_Standard.png',animal{kk}))
    saveas(gcf,sprintf('./saved_figs/%s_Standard.fig',animal{kk}))
end

f2 = figure;
subplot(1,2,1)
hold on
%plot3(bp_sine.param(:,1),bp_sine.param(:,2),gamma_norm_sine,'bo');
plot_surface_sine(gp_gamma_sine);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Gamma Power')
%view([-10,-10,8])
view(45,15)
subplot(1,2,2)
hold on
%plot3(bp_sine.param(:,1),bp_sine.param(:,2),theta_norm_sine,'bo');
plot_surface_sine(gp_theta_sine);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Theta Power')
suptitle('Sinusoid')
%view([-10,-10,8])
view(45,15)
f2.Position = [100 100 1000 600];

if save == 1
saveas(gcf,sprintf('./saved_figs/%s_Sine.png',animal{kk}))
saveas(gcf,sprintf('./saved_figs/%s_Sine.fig',animal{kk}))
end

f3 = figure;
subplot(1,2,1)
hold on
%plot3(bp_poisson.param(:,1),bp_poisson.param(:,2),gamma_norm_psn,'bo');
plot_surface(gp_gamma_psn);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Gamma Power')
view(45,15)

subplot(1,2,2)
hold on
%plot3(bp_poisson.param(:,1),bp_poisson.param(:,2),theta_norm_psn,'bo');
plot_surface(gp_theta_psn);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Theta Power')
view(45,15)
suptitle('Poisson Pulse')
f3.Position = [100 100 1000 600];
saveas(gcf,sprintf('./saved_figs/%s_Poisson.png',animal{kk}))
saveas(gcf,sprintf('./saved_figs/%s_Poisson.fig',animal{kk}))

train_width = [30 60];

f4 = figure;
subplot(2,2,1)
plot_gp_cross_sec(gp_gamma_npt,[2 3],[max(gp_gamma_npt.x_data(:,1)) train_width(1)/1000 .005])
xlabel('Train frequency (Hz)')
ylabel('Pulse frequency (Hz)')
title(sprintf('Train Width: %s ms',num2str(train_width(1))));
zlabel('Gamma power')

subplot(2,2,2)
plot_gp_cross_sec(gp_theta_npt,[2 3],[max(gp_gamma_npt.x_data(:,1)) train_width(1)/1000 .005])
xlabel('Train frequency (Hz)')
ylabel('Pulse frequency (Hz)')
zlabel('Theta power')
title(sprintf('Train Width: %s ms',num2str(train_width(1))));
subplot(2,2,3)
plot_gp_cross_sec(gp_gamma_npt,[2 3],[max(gp_gamma_npt.x_data(:,1)) train_width(2)/1000 .005])
xlabel('Train frequency (Hz)')
ylabel('Pulse frequency (Hz)')
zlabel('Gamma power')
title(sprintf('Train Width: %s ms',num2str(train_width(2))));
subplot(2,2,4)
plot_gp_cross_sec(gp_theta_npt,[2 3],[max(gp_gamma_npt.x_data(:,1)) train_width(2)/1000 .005])
xlabel('Train frequency (Hz)')
ylabel('Pulse frequency (Hz)')
zlabel('Theta power')
title(sprintf('Train Width: %s ms',num2str(train_width(2))));
suptitle('Nested Pulse Train')

f4.Position = [100 100 1000 1000];
if save == 1
saveas(gcf,sprintf('./saved_figs/%s_NPT.png',animal{kk}))
saveas(gcf,sprintf('./saved_figs/%s_NPT.fig',animal{kk}))
end



close all
% figure
% subplot(1,2,1)
% hold on
% histogram(theta_norm_std,'Normalization','pdf')
% %histogram(theta_norm_sine)
% histogram(theta_norm_npt,'Normalization','pdf')
% histogram(theta_norm_sine,'Normalization','pdf')
% legend({'Standard','Nested Pulse','Sine'})
% title('Theta Increase')
% 
% subplot(1,2,2)
% hold on
% histogram(gamma_norm_std,'Normalization','pdf')
% %histogram(theta_norm_sine)
% histogram(gamma_norm_npt,'Normalization','pdf')
% histogram(gamma_norm_sine,'Normalization','pdf')
% legend({'Standard','Nested Pulse','Sine'})
% title('Gamma Increase')


% 
% subplot(2,2,3)
% plot3(bandpower_data.param(:,1),bandpower_data.param(:,2),awake_gamma,'bo');
% xlabel('Amplitude')
% ylabel('Frequency')
% zlabel('Theta Power')
% subplot(2,2,4)
% plot_surface(gp_theta);
% xlabel('Amplitude')
% ylabel('Frequency')
% zlabel('Theta Power')


%difference in driving oscillations between standard pulse and
    %sinusoidal stimulation
    
    
end

%add Mark's repositories for GP modeling helper functions

%%%awake plot section


function plot_surface(gp_obj)
    t1 = linspace(gp_obj.lower_bound(1),gp_obj.upper_bound(1),50);
    t2 = linspace(gp_obj.lower_bound(2),gp_obj.upper_bound(2),50);
    t  = combvec(t1,t2)';
    state = repmat(.005,size(t,1),1); %plot

    [y_prediction, y_standard_deviation, fmu, fs2] = gp_obj.predict([t state]);
    y_pred_surf = reshape(y_prediction', [50,50]);
    surf(t1,t2,y_pred_surf');

end

function plot_surface_sine(gp_obj)
    t1 = linspace(gp_obj.lower_bound(1),gp_obj.upper_bound(1),50);
    t2 = linspace(gp_obj.lower_bound(2),gp_obj.upper_bound(2),50);
    t  = combvec(t1,t2)';
    %state = repmat(.005,size(t,1),1); %plot

    [y_prediction, y_standard_deviation, fmu, fs2] = gp_obj.predict([t]);
    y_pred_surf = reshape(y_prediction', [50,50]);
    surf(t1,t2,y_pred_surf');

end

function plot_gp_cross_sec(model,p_dim,p_vals)
    %p_dim: parameter space dimensions to plot on the cross section for NPT; 
    % p_dim = [2,3] means plot over train, pulse frequency
    % p_vals = [50, 40, 5] (values to keep other parameters constant - amplitude, train width, pulse width)
    %
    lower_bound = model.lower_bound;
    upper_bound = model.upper_bound;
    
    x1 = linspace(lower_bound(p_dim(1)),upper_bound(p_dim(1)),25);
    x2 = linspace(lower_bound(p_dim(2)),upper_bound(p_dim(2)),25);
    
    %amptemp = p_vals(1)*calib(1) + calib(2);
    amptemp = p_vals(1);
    
    input_space = combvec(x1,x2)';
    input_space = [amptemp*ones(size(input_space,1),1), input_space, p_vals(2)*ones(size(input_space,1),1) p_vals(3)*ones(size(input_space,1),1)];
    
    output = model.predict(input_space);
    to_plot = reshape(output,25,25);
    
    surf(x1,x2,to_plot)

end

function bp_norm = normalize_bp(bp_struct,band)

    switch band
        case 'gamma'
            bp_vals = bp_struct.lowgamma;
        case 'theta'
            bp_vals = bp_struct.theta;
    end
    
    bp_vals = mean(bp_vals,3);
    
    stim_times = find(bp_struct.stim_order == 1);
    
    if length(stim_times) == 1
        norm_temp = mean(bp_vals(:,1));
        bp_norm = (bp_vals(:,2)-norm_temp)/norm_temp*100;
        return
    end
    
    bp_norm = zeros(size(bp_vals,1),1);
    tempval = 2;
    norm_temp = mean(bp_vals(stim_times(1):stim_times(2),1));
    for kk = 1:size(bp_vals,1)
        if kk == stim_times(tempval) && tempval<length(stim_times)
            norm_temp = mean(bp_vals(stim_times(tempval):stim_times(tempval+1),1));
            tempval = tempval+1;
        else
            norm_temp = mean(bp_vals(stim_times(tempval):end,1));
        end
        bp_norm(kk) = (bp_vals(kk,2)-norm_temp)/norm_temp*100;
    end
    
end