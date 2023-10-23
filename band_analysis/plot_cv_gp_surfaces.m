
addpath('/Users/ERCOLE/Documents/Research/Data/Opto_GPs_2');
addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/tools/gpml-matlab-v3.6-2015-07-07'));

%addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/modeling'));
%addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/optogenetic_optimization_project'));


subject_ID = 'STV009';

switch subject_ID
    %get calibration, from manual input (stored in surgery log sheets)
    case 'ARN088'
        %opto_calibration = polyfit([3.8 3.98],[10 50]*(pi*.1^2),1);
        opto_calibration = polyfit([3.8 3.98],[10 50],1);
    case 'STV003'
        %opto_calibration = polyfit([.4 1.45],[10 50]*(pi*.1^2),1);
        opto_calibration = polyfit([.4 1.45],[10 50],1);
    case 'STV008'
        %opto_calibration = polyfit([.975 2.4375],[50 125]*(pi*.1^2),1);
        opto_calibration = polyfit([.975 2.4375],[50 125],1);
    case 'STV009'
        %opto_calibration = polyfit([.4 1.4],[10 50]*(pi*.1^2),1);
        opto_calibration = polyfit([.4 1.4],[10 50],1);
end

load(sprintf('%s_Awake_GP.mat',subject_ID));

gamma_std = gp.Standard.BandPower.Gamma;
theta_std = gp.Standard.BandPower.Theta;

gamma_psn = gp.Poisson.BandPower.Beta;
theta_psn = gp.Poisson.BandPower.Theta;

gamma_sine = gp.Sine.BandPower.Gamma;
theta_sine = gp.Sine.BandPower.Theta;

gamma_npt = gp.NPT.BandPower.Gamma;
theta_npt = gp.NPT.BandPower.Theta;

f1= figure(1);
subplot(1,2,1)
plot_surface_optoc(gamma_std,opto_calibration)
hold on
plot3((gamma_std.x_data(:,1)*opto_calibration(1)+opto_calibration(2)),gamma_std.x_data(:,2),(gamma_std.y_data-1)*100,'bo');
xlabel('Amplitude (mW/mm^2)')
ylabel('Frequency (Hz)')
zlabel('Gamma Power (% Change)')
title('Standard Pulse: Gamma')
set(gca,'FontSize',12)
view(45,15)

subplot(1,2,2)

plot_surface_optoc(theta_std,opto_calibration);
hold on
plot3(theta_std.x_data(:,1)*opto_calibration(1)+opto_calibration(2),theta_std.x_data(:,2),(theta_std.y_data-1)*100,'bo');
xlabel('Amplitude (mW/mm^2)')
ylabel('Frequency (Hz)')
zlabel('Theta Power (% Change)')
set(gca,'FontSize',12)
title('Standard Pulse: Theta')
view(45,15)
f1.Position = [100 100 1000 600];
saveas(gcf,sprintf('./saved_figs/%s_Standard.png',subject_ID))
saveas(gcf,sprintf('./saved_figs/%s_Standard.fig',subject_ID))

f2 = figure(2);
subplot(1,2,1)
hold on
plot3(gamma_psn.x_data(:,1)*opto_calibration(1)+opto_calibration(2),gamma_psn.x_data(:,2),(gamma_psn.y_data-1)*100,'bo');
plot_surface_optoc(gamma_psn,opto_calibration)
xlabel('Amplitude (mW/mm^2)')
ylabel('Frequency (Hz)')
zlabel('Gamma Power (% Change)')
set(gca,'FontSize',12)
title('Poisson Pulse: Gamma')
view(45,15)
subplot(1,2,2)
hold on
plot3(theta_psn.x_data(:,1)*opto_calibration(1)+opto_calibration(2),theta_psn.x_data(:,2),(theta_psn.y_data-1)*100,'bo');
plot_surface_optoc(theta_psn,opto_calibration);
xlabel('Amplitude (mW/mm^2)')
ylabel('Frequency (Hz)')
zlabel('Theta Power (% Change)')
set(gca,'FontSize',12)
title('Poisson Pulse: Theta')
view(45,15)
f2.Position = [100 100 1000 600];
saveas(gcf,sprintf('./saved_figs/%s_Poisson.png',subject_ID))
saveas(gcf,sprintf('./saved_figs/%s_Poisson.fig',subject_ID))

f3 = figure(3);

subplot(1,2,1)
hold on
plot3(gamma_sine.x_data(:,1)*opto_calibration(1)+opto_calibration(2),gamma_sine.x_data(:,2),(gamma_sine.y_data-1)*100,'bo');
plot_surface_sine_optoc(gamma_sine,opto_calibration)
xlabel('Amplitude (mW/mm^2)')
ylabel('Frequency (Hz)')
zlabel('Gamma Power (% Change)')
set(gca,'FontSize',12)
title('Sinusoid: Gamma')
view(45,15)
subplot(1,2,2)
hold on
plot3(theta_sine.x_data(:,1)*opto_calibration(1)+opto_calibration(2),theta_sine.x_data(:,2),(theta_sine.y_data-1)*100,'bo');
plot_surface_sine_optoc(theta_sine,opto_calibration);
xlabel('Amplitude (mW/mm^2)')
ylabel('Frequency (Hz)')
zlabel('Theta Power (% Change)')
set(gca,'FontSize',12)
title('Sinusoid: Theta')
view(45,15)
f3.Position = [100 100 1000 600];
saveas(gcf,sprintf('./saved_figs/%s_Sine.png',subject_ID))
saveas(gcf,sprintf('./saved_figs/%s_Sine.fig',subject_ID))

train_width = [30,60];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% uncomment for NPT GP's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% subplot(2,2,1)
% plot_gp_cross_sec(gamma_npt,[2 3],[max(gamma_npt.x_data(:,1)) train_width(1)/1000 .005],opto_calibration)
% xlabel('Train frequency (Hz)')
% ylabel('Pulse frequency (Hz)')
% title(sprintf('Train Width: %s ms',num2str(train_width(1))));
% zlabel('Gamma power')
% 
% subplot(2,2,2)
% plot_gp_cross_sec(theta_npt,[2 3],[max(gamma_npt.x_data(:,1)) train_width(1)/1000 .005],opto_calibration)
% xlabel('Train frequency (Hz)')
% ylabel('Pulse frequency (Hz)')
% zlabel('Theta power')
% title(sprintf('Train Width: %s ms',num2str(train_width(1))));
% subplot(2,2,3)
% plot_gp_cross_sec(gamma_npt,[2 3],[max(gamma_npt.x_data(:,1)) train_width(2)/1000 .005],opto_calibration)
% xlabel('Train frequency (Hz)')
% ylabel('Pulse frequency (Hz)')
% zlabel('Gamma power')
% title(sprintf('Train Width: %s ms',num2str(train_width(2))));
% subplot(2,2,4)
% plot_gp_cross_sec(theta_npt,[2 3],[max(gamma_npt.x_data(:,1)) train_width(2)/1000 .005],opto_calibration)
% xlabel('Train frequency (Hz)')
% ylabel('Pulse frequency (Hz)')
% zlabel('Theta power')
% title(sprintf('Train Width: %s ms',num2str(train_width(2))));
% suptitle('Nested Pulse Train')

%saveas(gcf,sprintf('/saved_figs/%s_NPT.png',subject_ID))
%saveas(gcf,sprintf('/saved_figs/%s_NPT.fig',subject_ID))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

function plot_surface_optoc(gp_obj,calib)
    t1 = linspace(gp_obj.lower_bound(1),gp_obj.upper_bound(1),50);
    t2 = linspace(gp_obj.lower_bound(2),gp_obj.upper_bound(2),50);
    t  = combvec(t1,t2)';
    state = repmat(.005,size(t,1),1); %plot

    [y_prediction, y_standard_deviation, fmu, fs2] = gp_obj.predict([t state]);
    y_pred_surf = reshape(y_prediction', [50,50]);
    surf(t1*calib(1)+calib(2),t2,100*(y_pred_surf-1)');

end

function plot_surface_sine_optoc(gp_obj,calib)
    t1 = linspace(gp_obj.lower_bound(1),gp_obj.upper_bound(1),50);
    t2 = linspace(gp_obj.lower_bound(2),gp_obj.upper_bound(2),50);
    t  = combvec(t1,t2)';
    %state = repmat(.005,size(t,1),1); %plot

    [y_prediction, y_standard_deviation, fmu, fs2] = gp_obj.predict([t]);
    y_pred_surf = reshape(y_prediction', [50,50]);
    surf(t1*calib(1)+calib(2),t2,100*(y_pred_surf-1)');

end

function plot_gp_cross_sec(model,p_dim,p_vals,calib)
    %p_dim: parameter space dimensions to plot on the cross section for NPT; 
    % p_dim = [2,3] means plot over train, pulse frequency
    % p_vals = [50, 40, 5] (values to keep other parameters constant - amplitude, train width, pulse width)
    %
    lower_bound = model.lower_bound;
    upper_bound = model.upper_bound;
    
    x1 = linspace(lower_bound(p_dim(1)),upper_bound(p_dim(1)),25);
    x2 = linspace(lower_bound(p_dim(2)),upper_bound(p_dim(2)),25);
    
    amptemp = p_vals(1)*calib(1) + calib(2);
    
    input_space = combvec(x1,x2)';
    input_space = [amptemp*ones(size(input_space,1),1), input_space, p_vals(2)*ones(size(input_space,1),1) p_vals(3)*ones(size(input_space,1),1)];
    
    output = model.predict(input_space);
    to_plot = reshape(output,25,25);
    
    surf(x1,x2,to_plot)

end