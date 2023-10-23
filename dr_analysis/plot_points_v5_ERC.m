%% Setup
clear
close all;
addpath('/Users/ERCOLE/Documents/Research/Repos/param_spaces_analysis/Violinplot-Matlab')
addpath('/Users/ERCOLE/Documents/Research/Repos/param_spaces_analysis')
version_string = '_v3_nptadd_30I';
% version_string = '_v2';
    % = '': Uses the original files that David used
    % = '_v2': Uses updated versions that I just made which include
        % recording data for STV003
    % = '_v2_30I': Version of files that exclude low-ampl stim parmeters
    %'_v3_nptadd_30I': best version

rescale_param = true;
rescale_f = 10;
load(sprintf('animals%s.mat',version_string))
animal = {'ARN088','STV009','STV003','STV008'};
kk = 1; % ARN088
% load(sprintf('ROI_%s%s.mat',animal{kk},version_string))

titles = {'UMAP','t-SNE','PCA'};
ii = 1; % Which DR algo to show: 1 = UMAP, 2 = tsne, 3 = PCA

% query_pt =  [.3,.50];%recording2  [-0.01 -0.18];%recording1  [-.78,.16];%branch pt   [-0.5 -0.1];%poisson pt  
%query_pt = [0.3,0.5;-0.01,-0.18;-0.78,0.16;-0.5,-0.1];
    %ORIGINAL LIST
    
%query_pt = [-.5,0;-0.25,-0.3;0.15,0.3;-0.05,0.05];
%query_pt = [0.3,0.5];
q_pt_to_show = 4;   %which point to gen figures to characterize parameters
   
    %FINALIZED QUERY POINT LIST FOR FIRST VERSION OF DR STRUCT
%query_pt = [-0.01,-0.18;0.3,0.5;-0.5,-0.1;-0.25,-0.3];

    %ADJUSTED QUERY POINT LIST AFTER RE-RUNNING DR PROCESSING
query_pt = [-0.02,0.014;0.30,0.45;-0.5,-0.07;-0.11,-0.48];   
spaces_to_plot = {[1,2],[1,2],[2,3],[2,3]};
% %Good query point list (might be out of date with newer file versions): 
%     %Arn88/poisson/tsne: [.55,.13]
%     %Arn88/recording1/tsne: [.55,.13]
%     %Arn88/recording2/tsne: [.01,-.16]
%     %Arn88/RecDoubleSine/tsne: [-.76,-.26]
%     %Arn88/NPT/tsne: [-.25,-.19]
%     %Arn88/Recording1/UMAP: [-.01,-.18]
%     %Arn88/StandardMix/UMAP: [-.78,.16]
% %     %Arn88/DoubleSine/UMAP: [.25,.59]
%     %STV009/Standard/UMAP: [.2,.5]
%     %STV009/Recording1/UMAP: [.11,-.76]
%     %STV009/Recording/t-SNE: [.36,-.37]
%     %STV009/Recording/tsne: [-.26,-.76]
%     %STV009/NPT/tsne: [-.62,-.29]
    
    
save_string = 'NPT';   %Name used when saving fig files
save = 1;   %binary flag to decide whether to save figs
    
Npts = 50;

spaces = {'Standard','Sine','Double Sine','Nested Pulse','Poisson','Recording'};

scatter_type = 'log';
colors = {[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.4660, 0.6740, 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],'k'};
%colors = {[0 0.30 0.7810],[0.8500 0.5 0.0980],[0.4940 0.1840 0.5560],[0.4660 0.9 0.1880],[0.9 0.0780 0.1840],'k'};
nic = [0.6350, 0.0780, 0.1840];
ic = [0, 0.5, 0];
sz = 5;

npts = 200;   %params for defining the grid for boundaries and stuff
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

[distances,idx] = sort(sqrt((data(:,1)-query_pt(q_pt_to_show,1)).^2+(data(:,2)-query_pt(q_pt_to_show,2)).^2));

sorted_params = param_stack(idx(1:Npts));
sorted_types = param_types(idx(1:Npts));

% scatter(data(idx(1:Npts),1),data(idx(1:Npts),2),50,[0, 0.5, 0],'filled')
scatter(query_pt(:,1),query_pt(:,2),140,'o','MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',4)
hold off
%title(titles{ii},'FontSize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',16)
legend({'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior','Query Point'},'Location','Southeast');
xlabel('Dim. 1')
ylabel('Dim. 2')

if save == 1
    %saveas(gcf,sprintf('Figures/%s_%s_scatter_%.2f_%.2f.png',animal{kk},titles{ii},query_pt(1),query_pt(2)))
    %saveas(gcf,sprintf('./query_figs/%s_%s_pt%d_Scatter2.png',animal{kk},titles{ii},q_pt_to_show))
    exportgraphics(gcf,'UMAP_w_queries.png','Resolution',300)
    %saveas(gcf,sprintf('./query_figs/%s_%s_pt%d_Scatter2.svg',animal{kk},titles{ii},q_pt_to_show))
end

%figure
figure('Position',[100,100,300,700])
subplot(2,2,[1 2])
lim = 0;
cutoff_freq = 5;
space = {};
labels = {};
% subplot(1,2,1)
% 
% for j = 1:Npts
%     plot((cutoff_freq:size(psd,2))*0.3049,log10(psd(idx(j),cutoff_freq:end)),'Color',colors{grp(idx(j))},'LineWidth',1)
%     hold on
% %     if max(psd(idx(j),cutoff_freq:end)) > lim
% %         lim = max(psd(idx(j),cutoff_freq:end));
% %     end
%     space{j,1} = spaces{grp(idx(j))};
%     labels{j,1} = 'Frequency';
%     labels{j,2} = 'Amplitude';
%     labels{j,3} = '2nd Frequency';
% end
% hold off
% %ylim([-lim*.1,lim])
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)
set(gca,'fontsize',14)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (dB)','FontSize',14)

rec_inds = find(contains(param_types,'Recording'));

psd_mean = log10(mean(psd(idx(1:Npts),cutoff_freq:end)));
psd_serr = std(log10(psd(idx(1:Npts),cutoff_freq:end)))/sqrt(Npts);
fr_patch = [(cutoff_freq:size(psd,2))*0.3049, fliplr(cutoff_freq:size(psd,2))*0.3049];
psd_patch = [psd_mean-psd_serr, fliplr(psd_mean+psd_serr)];

psd_mean_all = log10(mean(psd(:,cutoff_freq:end)));
psd_serr_all = std(log10(psd(:,cutoff_freq:end)))/sqrt(Npts);
psd_patch_all = [psd_mean_all-psd_serr_all, fliplr(psd_mean_all+psd_serr_all)];

psd_mean_rec = log10(mean(psd(rec_inds,cutoff_freq:end)));
psd_serr_rec = std(log10(psd(rec_inds,cutoff_freq:end)))/sqrt(Npts);
psd_patch_rec = [psd_mean_rec-psd_serr_rec, fliplr(psd_mean_rec+psd_serr_rec)];

% subplot(1,2,2)
hold on
patch(fr_patch,psd_patch,[1 .4 .08]*.85);
%patch(fr_patch,psd_patch_all,[0.5,0.5,0.5],'FaceAlpha',0.3)   
    %^plot average of all data for comparison
patch(fr_patch,psd_patch_rec,[0.5,0.5,0.5],'FaceAlpha',0.3)  
legend({'Queried Points','Behavior'})
axis([0 50 -8.5 -3]);
set(gca,'fontsize',13)
%suptitle(sprintf('Power Spectral Densities around point (%.2f, %.2f)',query_pt(q_pt_to_show,1),query_pt(q_pt_to_show,2)))
if save == 1
    %saveas(gcf,sprintf('Figures/%s_%s_psd_%.2f_%.2f.png',animal{kk},titles{ii},query_pt(1),query_pt(2)))
%     saveas(gcf,sprintf('./query_figs/%s_%s_pt%d_PSD2.png',animal{kk},titles{ii},q_pt_to_show))
%     saveas(gcf,sprintf('./query_figs/%s_%s_pt%d_PSD2.fig',animal{kk},titles{ii},q_pt_to_show))
end

% figure()
sort_param_counts = tabulate(sorted_types);
% bar(cat(1,sort_param_counts{:,2}))
% set(gca,'XtickLabel',sort_param_counts(:,1))
% ylabel('Count')
% title('Parameter Space Count')

fprintf('Parameter Space Counts:\n')
for kprint = 1:size(sort_param_counts,1)
    fprintf('    %s: %d trials\n',sort_param_counts{kprint,1},sort_param_counts{kprint,2})
end

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
% if save == 1
%     saveas(gcf,sprintf('query_figures/%s_%s_%s_%d_ParamType.png',animal{kk},save_string,titles{ii},q_pt_to_show))
%     saveas(gcf,sprintf('query_figures/%s_%s_%s_%d_ParamType.fig',animal{kk},save_string,titles{ii},q_pt_to_show))
% end

%figure('Position', [100 100 1000 400])

loopind = 1;
for k2 = spaces_to_plot{q_pt_to_show}%1:max(sort_params_inds)
    
    %subplot(1,3,loopind)
    subplot(2,2,loopind+2)
    if strcmp(sort_params_unique{k2},'Recording')
       continue 
    elseif strcmp(sort_params_unique{k2},'Standard') || strcmp(sort_params_unique{k2},'Poisson')
        
        labels = {'Ampl.','Freq.','Pulse Width'};
        temp = sort_params_vals{k2};
        temp(:,3) = temp(:,3)*1000;
        %v = violinplot(temp,labels);
        data_plot = temp;
        ylim([0 50])
        if strcmp(sort_params_unique{k2},'Standard')
            coltemp = colors{1};
        end
        if strcmp(sort_params_unique{k2},'Poisson')
            coltemp = colors{5};
        end
    elseif strcmp(sort_params_unique{k2},'Nested Pulse')
        labels = {'Ampl.','Train Freq.','Pulse Freq.','Train Width','Pulse Width'};
        temp = sort_params_vals{k2};
        temp(:,4:5) = temp(:,4:5)*1000;
        if rescale_param
           temp(:,2) = temp(:,2)*rescale_f;
           labels{2} = 'Train Freq. (x10)';
           temp(:,5) = temp(:,5)*rescale_f;
           labels{5} = 'Pulse Width (x10)';
        end
        data_plot = temp;
        coltemp = colors{4};
        for k3 = 1:size(data_plot,2)
            Violin(data_plot(:,k3),k3,'ViolinColor',coltemp)
        end
        xticks(1:size(data_plot,2))
        xticklabels(labels)
        %v = violinplot(temp,labels);
        ylim([0 100])
    elseif strcmp(sort_params_unique{k2},'Sine')
        coltemp = colors{2};
        labels = {'Ampl.','Freq.'};
        data_plot = sort_params_vals{k2};
        %v = violinplot(sort_params_vals{k2},labels);
        ylim([0 50])
    elseif strcmp(sort_params_unique{k2},'Double Sine')
        coltemp = colors{3};
        labels = {'Ampl.','Ratio','Main Freq.','2nd Freq.'};
        data_plot = sort_params_vals{k2};
        for km = 1:size(data_plot,1)
            if (data_plot(km,2)) < 0.5   %Switch order of double sine parameters so "main"frequency
                                    % is ordered first
                temp = data_plot(km,3);
                data_plot(km,3) = data_plot(km,4);
                data_plot(km,4) = temp;    
            end
        end
        if rescale_param
           data_plot(:,2) =  data_plot(:,2)*10;
           labels{2} = 'Ratio (x10)';
        end
        
    end
    if ~strcmp(sort_params_unique{k2},'Nested Pulse')
        for k3 = 1:size(data_plot,2)
            Violin(data_plot(:,k3),k3,'ViolinColor',coltemp)
        end
        xticks(1:size(data_plot,2))
        xticklabels(labels)
        %v = violinplot(data_plot,labels);
        ylim([0 50])
    end
    title(sort_params_unique{k2})
    xtickangle(45)
    loopind = loopind+1;
    set(gca,'fontsize',13)
end

if save == 1
    %saveas(gcf,sprintf('./query_figs/%s_%s_%d_ParamVals2.png',animal{kk},titles{ii},q_pt_to_show))
    saveas(gcf,sprintf('./query_figs/%s_%s_%d_ParamVals2.svg',animal{kk},titles{ii},q_pt_to_show))
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
