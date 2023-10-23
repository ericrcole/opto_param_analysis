%% Setup
clear all;
addpath('Violinplot-Matlab')

version_string = '_v2_30I';
    % = '': Uses the original files that David used
    % = '_v2': Uses updated versions that I just made which include
        % recording data for STV003
    % = '_v2_30I': Version of files that exclude low-ampl stim parmeters

load(sprintf('animals%s.mat',version_string))
animal = {'ARN088','STV009','STV003','STV008'};
kk = 1; % ARN088
load(sprintf('ROI_%s%s.mat',animal{kk},version_string))

titles = {'UMAP','t-SNE','PCA'};
ii = 1; % Which DR algo to show: 1 = UMAP, 2 = tsne, 3 = PCA

query_pt = [-.62,-.29];%[-0.90, -0.07];
    %Good query point list (might be out of date with newer file versions): 
    %Arn88/poisson/tsne: [.55,.13]
    %Arn88/recording1/tsne: [.55,.13]
    %Arn88/recording2/tsne: [.01,-.16]
    %Arn88/RecDoubleSine/tsne: [-.76,-.26]
    %Arn88/NPT/tsne: [-.25,-.19]
    %Arn88/Recording1/UMAP: [-.01,-.18]
    %Arn88/StandardMix/UMAP: [-.78,.16]
    %Arn88/DoubleSine/UMAP: [.25,.59]
    %STV009/Standard/UMAP: [.2,.5]
    %STV009/Recording1/UMAP: [.11,-.76]
    %STV009/Recording/t-SNE: [.36,-.37]
    %STV009/Recording/tsne: [-.26,-.76]
    %STV009/NPT/tsne: [-.62,-.29]
    
    
save_string = 'NPT';   %Name used when saving fig files
save = 0;   %binary flag to decide whether to save figs
    
Npts = 50;

spaces = {'Standard','Sine','Double Sine','Nested Pulse','Poisson','Recording'};

scatter_type = 'log';
colors = {[0 0.30 0.7810],[0.8500 0.5 0.0980],[0.4940 0.1840 0.5560],[0.4660 0.9 0.1880],[0.9 0.0780 0.1840],'k'};
nic = [0.6350, 0.0780, 0.1840];
ic = [0, 0.5, 0];
sz = 5;

npts = 200;
kernel_width = 0.05; %0.03
prob_thresh = 0.01;


endloop = 6;

xmax = zeros(3,1);


param_stack = {}; 
param_types = {};
for k1 = 1:length(animals{kk}.clusters{1}.param)
    param_stack = [param_stack; num2cell(animals{kk}.clusters{1}.param{k1},2)];
    celltemp    = cell(size(animals{kk}.clusters{1}.param{k1},1),1);
    celltemp(:) = {spaces{k1}};
    param_types = [param_types; celltemp];
end
%animals{kk}.clusters{1}.param

xlims = zeros(1,2);
ylims = zeros(1,2);
%% Plot points
for kl = 1:endloop
    dat = animals{kk}.clusters{ii}.data{kl};
    xlims = [xlims; [min(dat(:,1)) max(dat(:,1))]];
    ylims = [ylims; [min(dat(:,2)) max(dat(:,2))]];
end
xmax(ii) = max(abs(xlims(:)));
ymax(ii) = max(abs(ylims(:)));
% Make Grid over which to evaluate the RBF kernel density
xlimits = [-1.1 1.1];
ylimits = [-1.1 1.1];
%xlimits = [floor(min(xlims(:,1)))-3 ceil(max(xlims(:,2)))+3];
%ylimits = [floor(min(ylims(:,1)))-3 ceil(max(ylims(:,2)))+3];
xpts = linspace(xlimits(1), xlimits(2),npts);
ypts = linspace(ylimits(1), ylimits(2),npts);
[X,Y] = meshgrid(xpts,ypts);
xi = [X(:), Y(:)];
XY{ii} = xi;
volume_ovl = [];

figure()
X = reshape(XY{ii}(:,1),npts,npts);
Y = reshape(XY{ii}(:,2),npts,npts);
data = [];
psd = [];
grp = [];
freqs = [];
freqs2 = [];
amps = [];
for kl = 1:endloop
    hold on
    dat = animals{kk}.clusters{ii}.data{kl}./[xmax(ii) ymax(ii)];
    power_sd = animals{kk}.clusters{ii}.psd{kl};
    group = ones(length(dat),1)*kl;
    param = animals{kk}.clusters{ii}.param{kl};

    f = mvksdensity(dat,xi,'Bandwidth',kernel_width);
    f_reshape = reshape(f,size(X));
    % surf(X,Y,f_reshape)
    idx = find(f_reshape > prob_thresh*max(f));
    idxs{kl,ii} = idx;

    scatter(dat(:,1),dat(:,2),sz,colors{kl},'filled')
    %h = plot(cloud{kl,ii}(rbf_boundary,1),cloud{kl,ii}(rbf_boundary,2),colors{kl},'LineWidth',2);
    %h = scatter(cloud(rbf_boundary,1),cloud(rbf_boundary,2),colors{kl},'LineWidth',2);
    %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %spaces = {'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'};
    switch kl
        case 1 % Standard
            freq = param(:,2);
            amp = param(:,1);
            pulse_width = param(:,3);
            freq2 = zeros(length(dat),1);
        case 2 % Sine
            freq = param(:,2);
            amp = param(:,1);
            freq2 = zeros(length(dat),1);
        case 3 % Sine_2N
            freq = param(:,3);
            amp = param(:,1);
            freq2 = param(:,4);
        case 4 % NPT
            freq = param(:,3);
            amp = param(:,1);
            freq2 = zeros(length(dat),1);
        case 5 % Poisson
            freq = param(:,2);
            amp = param(:,1);
            freq2 = zeros(length(dat),1);
        case 6 % Recording
            freq = zeros(length(dat),1);
            amp = zeros(length(dat),1);
            freq2 = zeros(length(dat),1);
    end
    data = [data; dat];
    psd = [psd; power_sd];
    grp = [grp; group];
    freqs = [freqs; freq];
    amps = [amps; amp];
    freqs2 = [freqs2; freq2];
end

[distances,idx] = sort(sqrt((data(:,1)-query_pt(1)).^2+(data(:,2)-query_pt(2)).^2));

sorted_params = param_stack(idx(1:Npts));
sorted_types = param_types(idx(1:Npts));

scatter(data(idx(1:Npts),1),data(idx(1:Npts),2),50,[0, 0.5, 0],'filled')
scatter(query_pt(1),query_pt(2),140,'ko','filled')
hold off
title(titles{ii},'FontSize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)
legend({'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior','Neighbors','Query Point'});
if save == 1
    %saveas(gcf,sprintf('Figures/%s_%s_scatter_%.2f_%.2f.png',animal{kk},titles{ii},query_pt(1),query_pt(2)))
    saveas(gcf,sprintf('query_figures/%s_%s_%s_Scatter.png',animal{kk},save_string,titles{ii}))
    saveas(gcf,sprintf('query_figures/%s_%s_%s_Scatter.fig',animal{kk},save_string,titles{ii}))
end

figure()
lim = 0;
cutoff_freq = 5;
space = {};
labels = {};
subplot(1,2,1)

for j = 1:Npts
    plot((cutoff_freq:size(psd,2))*0.3049,log10(psd(idx(j),cutoff_freq:end)),'Color',colors{grp(idx(j))},'LineWidth',1)
    hold on
%     if max(psd(idx(j),cutoff_freq:end)) > lim
%         lim = max(psd(idx(j),cutoff_freq:end));
%     end
    space{j,1} = spaces{grp(idx(j))};
    labels{j,1} = 'Frequency';
    labels{j,2} = 'Amplitude';
    labels{j,3} = '2nd Frequency';
end
hold off
%ylim([-lim*.1,lim])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)
xlabel('Frequency (Hz)','FontSize',16)

psd_mean = log10(mean(psd(idx(1:Npts),cutoff_freq:end)));
psd_serr = std(log10(psd(idx(1:Npts),cutoff_freq:end)))/sqrt(Npts);
fr_patch = [(cutoff_freq:size(psd,2))*0.3049, fliplr(cutoff_freq:size(psd,2))*0.3049];
psd_patch = [psd_mean-psd_serr, fliplr(psd_mean+psd_serr)];

subplot(1,2,2)
patch(fr_patch,psd_patch,[1 .4 .08]*.85);
suptitle(sprintf('Power Spectral Densities around point (%.2f, %.2f)',query_pt(1),query_pt(2)))
if save == 1
    %saveas(gcf,sprintf('Figures/%s_%s_psd_%.2f_%.2f.png',animal{kk},titles{ii},query_pt(1),query_pt(2)))
    saveas(gcf,sprintf('query_figures/%s_%s_%s_PSD.png',animal{kk},save_string,titles{ii}))
    saveas(gcf,sprintf('query_figures/%s_%s_%s_PSD.fig',animal{kk},save_string,titles{ii}))
end

figure()
sort_param_counts = tabulate(sorted_types)
bar(cat(1,sort_param_counts{:,2}))
set(gca,'XtickLabel',sort_param_counts(:,1))
ylabel('Count')
title('Parameter Space Count')

[sort_params_unique,~,sort_params_inds] = unique(sorted_types);
sort_params_vals = cell(size(sort_params_unique));
for k2 = 1:max(sort_params_inds)
    if strcmp(sort_params_unique{k2},'Recording')
       continue 
    end
    temp = idx(1:Npts);
    temp_log = temp(sort_params_inds == k2);
    sort_params_vals{k2} = [sort_params_vals{k2}; cat(1,param_stack{temp_log})];
    
end
if save == 1
    saveas(gcf,sprintf('query_figures/%s_%s_%s_ParamType.png',animal{kk},save_string,titles{ii}))
    saveas(gcf,sprintf('query_figures/%s_%s_%s_ParamType.fig',animal{kk},save_string,titles{ii}))
end

figure()
for k2 = 1:max(sort_params_inds)
    
    subplot(2,3,k2)
    if strcmp(sort_params_unique{k2},'Recording')
       continue 
    elseif strcmp(sort_params_unique{k2},'Standard') || strcmp(sort_params_unique{k2},'Poisson')
        temp = sort_params_vals{k2};
        temp(:,3) = temp(:,3)*1000;
        v = violinplot(temp,labels);
    elseif strcmp(sort_params_unique{k2},'NPT')
        temp = sort_params_vals{k2};
        temp(:,4:5) = temp(:,4:5)*1000;
        v = violinplot(temp,labels);
    else
        v = violinplot(sort_params_vals{k2},labels);
    end
    
    title(sort_params_unique{k2})
end
suptitle('Parameter values around query point')

if save == 1
    saveas(gcf,sprintf('query_figures/%s_%s_%s_ParamVals.png',animal{kk},save_string,titles{ii}))
    saveas(gcf,sprintf('query_figures/%s_%s_%s_ParamVals.fig',animal{kk},save_string,titles{ii}))
end
    %suptitle(sprintf('Parameters around point (%.2f, %.2f)',query_pt(1),query_pt(2)),'FontSize',16)

% figure()
% frequencies = freqs(idx(1:Npts));
% amplitudes = amps(idx(1:Npts));
% frequencies2 = freqs2(idx(1:Npts));
% %labels = cellstr(['Frequency','Amplitude']);
% grps = grp(idx(1:Npts));
% v = violinplot([frequencies,amplitudes,frequencies2],labels);
% title(sprintf('Parameters around point (%.2f, %.2f)',query_pt(1),query_pt(2)),'FontSize',16)
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)
% saveas(gcf,sprintf('Figures/%s_%s_params_%.2f_%.2f.png',animal{kk},titles{ii},query_pt(1),query_pt(2)))
% 
