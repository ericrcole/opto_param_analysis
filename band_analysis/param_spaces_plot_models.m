%model_NPT_data.m
%
%Trains multi-dimensional GP model for theta power vs. NPT parameters
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/tools/gpml-matlab-v3.6-2015-07-07'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/modeling'));
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/optogenetic_optimization_project'));

cd('/Users/ERCOLE/Documents/Research/Data');

load('ARN088_Bandpower_NPT.mat');

theta_NPT = mean(bandpower_data.theta(:,2),3);
gamma_NPT = mean(bandpower_data.lowgamma(:,2),3);
param = bandpower_data.param;

histogram(theta_NPT)

p_outliers = .04; %top percent of trials to discard as outliers
n_outliers = round(p_outliers*length(theta_NPT));
[~,out_inds] = maxk(theta_NPT,n_outliers);

theta_NPT(out_inds) = [];
param(out_inds,:) = [];
gamma_NPT(out_inds) = [];

param(:,4:5) = param(:,4:5)*1000;

amp_inds = and(param(:,1) > 30, param(:,1) < 50);

figure
subplot(1,2,1)
scatter3(param(amp_inds,2),param(amp_inds,3),theta_NPT(amp_inds));
xlabel('Train frequency')
ylabel('Pulse frequency')
title('Theta')
subplot(1,2,2)
scatter3(param(amp_inds,2),param(amp_inds,3),gamma_NPT(amp_inds));
xlabel('Train frequency')
ylabel('Pulse frequency')
title('Gamma')

gp_theta = gp_object();
gp_gamma = gp_object();

gp_theta.initialize_data(param,theta_NPT);
gp_gamma.initialize_data(param,gamma_NPT);

% gp_theta.minimize(10);
% gp_gamma.minimize(10);

figure
suptitle('Nested Pulse Train Stim')
subplot(2,2,1)
plot_gp_cross_sec(gp_theta,[2 3],[40 40 5]);
xlabel('Train frequency')
ylabel('Pulse frequency')
title('Theta')
subplot(2,2,2)
plot_gp_cross_sec(gp_gamma,[2 3],[40 40 5]);
xlabel('Train frequency')
ylabel('Pulse frequency')
title('Gamma')
subplot(2,2,3)
plot_gp_cross_sec(gp_theta,[2 3],[40 70 5]);
xlabel('Train frequency')
ylabel('Pulse frequency')
title('Theta')
subplot(2,2,4)
plot_gp_cross_sec(gp_gamma,[2 3],[40 70 5]);
xlabel('Train frequency')
ylabel('Pulse frequency')
title('Gamma')

clear bandpower_data

load('ARN088_Bandpower_Sine_2N.mat');

theta_sine = mean(bandpower_data.theta(:,2),3);
gamma_sine = mean(bandpower_data.lowgamma(:,2),3);
param = bandpower_data.param;

gp_theta_sine = gp_object();
gp_gamma_sine = gp_object();

gp_theta_sine.initialize_data(param,theta_sine);
gp_gamma_sine.initialize_data(param,gamma_sine);

figure
suptitle('Summed Sinusoids')
subplot(2,2,1)
plot_gp_cross_sine(gp_theta_sine,[3 4],[30 .33]);
xlabel('Frequency 1')
ylabel('Frequency 2')
title('Theta')
subplot(2,2,2)
plot_gp_cross_sine(gp_gamma_sine,[3 4],[30 .66]);
xlabel('Frequency 1')
ylabel('Frequency 2')
title('Gamma')
subplot(2,2,3)
plot_gp_cross_sine(gp_theta_sine,[3 4],[30 .33]);
xlabel('Frequency 1')
ylabel('Frequency 2')
title('Theta')
subplot(2,2,4)
plot_gp_cross_sine(gp_gamma_sine,[3 4],[30 .66]);
xlabel('Frequency 1')
ylabel('Frequency 2')
title('Gamma')

load('ARN088_Bandpower_Poisson.mat');

theta_psn = mean(bandpower_data.theta(:,2),3);
gamma_psn = mean(bandpower_data.lowgamma(:,2),3);
param = bandpower_data.param;

gp_theta_psn = gp_object();
gp_gamma_psn = gp_object();

gp_theta_psn.initialize_data(param,theta_psn);
gp_gamma_psn.initialize_data(param,gamma_psn);

figure
suptitle('Poisson Stim')
subplot(2,2,1)
plot3(bandpower_data.param(:,1),bandpower_data.param(:,2),gamma_psn,'bo');
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Gamma Power')
title('Gamma')
subplot(2,2,2)
plot_surface(gp_gamma_psn);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Gamma Power')

subplot(2,2,3)
plot3(bandpower_data.param(:,1),bandpower_data.param(:,2),theta_psn,'bo');
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Theta Power')
title('Theta')
subplot(2,2,4)
plot_surface(gp_theta_psn);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Theta Power')

load('ARN088_Bandpower_Standard.mat');

theta_std = mean(bandpower_data.theta(:,2),3);
gamma_std = mean(bandpower_data.lowgamma(:,2),3);
param = bandpower_data.param;

n_outliers = round(p_outliers*length(theta_std));
[~,out_inds] = maxk(theta_std,n_outliers);

theta_std(out_inds) = [];
param(out_inds,:) = [];
gamma_std(out_inds) = [];

gp_theta_std = gp_object();
gp_gamma_std = gp_object();

gp_theta_std.initialize_data(param,theta_std);
gp_gamma_std.initialize_data(param,gamma_std);

figure
suptitle('Standard Pulse Stim')
subplot(2,2,1)
plot3(param(:,1),param(:,2),gamma_std,'bo');
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Gamma Power')
title('Gamma Power')
subplot(2,2,2)
plot_surface(gp_gamma_std);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Gamma')

subplot(2,2,3)
plot3(param(:,1),param(:,2),theta_std,'bo');
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Theta Power')
title('Theta')
subplot(2,2,4)
plot_surface(gp_theta_std);
xlabel('Amplitude')
ylabel('Frequency')
zlabel('Theta Power')

function plot_gp_cross_sec(model,p_dim,p_vals)
    lower_bound = model.lower_bound;
    upper_bound = model.upper_bound;
    
    x1 = linspace(lower_bound(p_dim(1)),upper_bound(p_dim(1)),25);
    x2 = linspace(lower_bound(p_dim(2)),upper_bound(p_dim(2)),25);
    input_space = combvec(x1,x2)';
    input_space = [p_vals(1)*ones(size(input_space,1),1), input_space, p_vals(2)*ones(size(input_space,1),1) p_vals(3)*ones(size(input_space,1),1)];
    
    output = model.predict(input_space);
    to_plot = reshape(output,25,25);
    
    surf(x1,x2,to_plot)

end

function plot_gp_cross_sine(model,p_dim,p_vals)
    lower_bound = model.lower_bound;
    upper_bound = model.upper_bound;
    
    x1 = linspace(lower_bound(p_dim(1)),upper_bound(p_dim(1)),25);
    x2 = linspace(lower_bound(p_dim(2)),upper_bound(p_dim(2)),25);
    input_space = combvec(x1,x2)';
    input_space = [p_vals(1)*ones(size(input_space,1),1), p_vals(2)*ones(size(input_space,1),1) input_space];
    
    output = model.predict(input_space);
    to_plot = reshape(output,25,25);
    
    surf(x1,x2,to_plot)

end

function plot_surface(gp_obj)
    t1 = linspace(gp_obj.lower_bound(1),gp_obj.upper_bound(1),50);
    t2 = linspace(gp_obj.lower_bound(2),gp_obj.upper_bound(2),50);
    t  = combvec(t1,t2)';
    state = repmat(.005,size(t,1),1); %plot

    [y_prediction, y_standard_deviation, fmu, fs2] = gp_obj.predict([t state]);
    y_pred_surf = reshape(y_prediction', [50,50]);
    surf(t1,t2,y_pred_surf');

end
