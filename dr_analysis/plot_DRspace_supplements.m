


version_string = '_v2_30I';
load(sprintf('animals%s.mat',version_string))
addpath('Violinplot-Matlab')
animal = {'ARN088','STV009','STV003','STV008'};
ii_def = 1; % Use UMAP
plot_animal = 1; % only plot ARN088
titles = {'UMAP','t-SNE','PCA'};

spaces = {'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'};

addpath('/Volumes/OPTO_PROC/OptoStim')

scatter_type = 'log';
%colors = {'g','m','c','r','b','k'};
colors = {[0 0.30 0.7810],[0.8500 0.5 0.0980],[0.4940 0.1840 0.5560],[0.4660 0.9 0.1880],[0.9 0.0780 0.1840],'k'};
nic = [0.6350, 0.0780, 0.1840];
ic = [0, 0.5, 0];
sz = 10;

npts = 200;
kernel_width = 0.05; %0.03
prob_thresh = 0.01;


xmax = zeros(3,1);

f1 = figure('Position',  [100, 100, 800, 800]);
%showing the DR plots for all subjects, DR methods in one figure

% opto_slopes = zeros(size(animal));
% opto_ints = zeros(size(animal));

mean_gamma = zeros(1,4);
mean_gamma_low = zeros(1,4);
separability = zeros(3,4);

for kk = 1:4
    for ii = 1:3
        load(sprintf('%s_BandpowerPhase_Awake_Standard.mat',animal{kk}));
        amp_max = max(bandpower_data.param(:,1));
        amp_min = min(bandpower_data.param(:,1));
        param_inds_low = find(and(bandpower_data.param(:,1) == amp_min,bandpower_data.param(:,2) == 35));
        param_inds = find(and(bandpower_data.param(:,1) == amp_max,bandpower_data.param(:,2) == 35));
        gamma_vals = mean(bandpower_data.lowgamma,3);
        mean_gamma(kk) = mean((gamma_vals(param_inds,2) - gamma_vals(param_inds,1))./gamma_vals(param_inds,1)*100);
        mean_gamma_low(kk) = mean((gamma_vals(param_inds_low,2) - gamma_vals(param_inds_low,1))./gamma_vals(param_inds_low,1)*100);
        
        %     switch animal{kk}
        %         case 'ARN088'
        %             opto_calibration = polyfit([3.8 3.98],[10 50]*(pi*.1^2),1);
        %         case 'STV003'
        %             opto_calibration = polyfit([.4 1.45],[10 50]*(pi*.1^2),1);
        %         case 'STV008'
        %             opto_calibration = polyfit([.975 2.4375],[50 125]*(pi*.1^2),1);
        %         case 'STV009'
        %             opto_calibration = polyfit([.4 1.4],[10 50]*(pi*.1^2),1);
        %     end
        
        endloop = 6;
        
        xlims = zeros(1,2);
        ylims = zeros(1,2);
        
        % kl is clusters
        for kl = 1:endloop
            dat = animals{kk}.clusters{ii}.data{kl};
            xlims = [xlims; [min(dat(:,1)) max(dat(:,1))]];
            ylims = [ylims; [min(dat(:,2)) max(dat(:,2))]];
        end
        xmax(kk) = max(abs(xlims(:)));
        ymax(kk) = max(abs(ylims(:)));
        
        % Make Grid over which to evaluate the RBF kernel density
        xlimits = [-1.1 1.1];
        ylimits = [-1.1 1.1];
        %xlimits = [floor(min(xlims(:,1)))-3 ceil(max(xlims(:,2)))+3];
        %ylimits = [floor(min(ylims(:,1)))-3 ceil(max(ylims(:,2)))+3];
        xpts = linspace(xlimits(1), xlimits(2),npts);
        ypts = linspace(ylimits(1), ylimits(2),npts);
        [X,Y] = meshgrid(xpts,ypts);
        xi = [X(:), Y(:)];
        XY{kk} = xi;
        
        subplot(4,3,3*kk+ii-3)
        for kl = 1:endloop
            
            hold on
            
            dat = animals{kk}.clusters{ii}.data{kl}./[xmax(kk) ymax(kk)];
            %plot_error_ellipse(animals{kk}.clusters{ii}.dat{kl},colors{kl},false)
            
            f = mvksdensity(dat,xi,'Bandwidth',kernel_width);
            f_reshape = reshape(f,size(X));
            % surf(X,Y,f_reshape)
            idx = find(f_reshape > prob_thresh*max(f));
            idxs{kl,kk} = idx;
            
            scatter(dat(:,1),dat(:,2),sz,colors{kl},'filled')
            
            %h = plot(cloud{kl,ii}(rbf_boundary,1),cloud{kl,ii}(rbf_boundary,2),colors{kl},'LineWidth',2);
            %h = scatter(cloud(rbf_boundary,1),cloud(rbf_boundary,2),colors{kl},'LineWidth',2);
            %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            hold off
            ROI{kl,kk} = zeros(size(f_reshape));
            ROI{kl,kk}(idx) = 1;
            
        end
        
        if rem(3*kk+ii-3,3) == 1
            ylabel(sprintf('Subject %d',kk),'FontSize',14)
        end
        if 3*kk+ii-3<4
            title(titles{ii},'FontSize',14)
        end
        %title(sprintf('Subject %d',kk),'FontSize',16)
        %xlabel('Dimension 1','FontSize',16)
        %ylabel('Dimension 2','FontSize',16)
        
        %     legend({'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'},'Location','Northeast');
        %     else
        %     legend({'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'},'Location','Southeast');
        %     end
        %     xlim(xlimits)
        %     ylim(ylimits)
        
        %exportgraphics(f,sprintf('Figures/%s_data.fig',animal{kk}),'Resolution',300)
        
        %     if kk == plot_animal
        %         subplot(2,3,2)
        %         for kl = 1:endloop
        %             hold on
        %             dat = animals{kk}.clusters{ii}.data{kl}./[xmax(kk) ymax(kk)];
        %             img = zeros(size(X));
        %             img(idxs{kl,kk}) = 1;
        %             bounds = bwboundaries(img);
        %             for k = 1 : length(bounds)
        %                 thisBoundary = bounds{k};
        %                 x = thisBoundary(:, 2);
        %                 y = thisBoundary(:, 1);
        %                 vector = [];
        %                 for l = 1:length(x)
        %                     vector = [vector; [X(x(l),y(l)) Y(x(l),y(l))]];
        %                 end
        %                 h = plot(vector(:,2),vector(:,1),'Color',colors{kl},'Linewidth',2);
        %                 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %             end
        %             scatter(dat(:,1),dat(:,2),sz,colors{kl},'filled')
        %         end
        %         legend({'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'},'Location','Southeast');
        %         xlabel('Dimension 1','FontSize',16)
        %         ylabel('Dimension 2','FontSize',16)
        %         title(sprintf('Boundaries of Parameter Spaces'),'FontSize',16)
        %         xlim(xlimits)
        %         ylim(ylimits)
        %         hold off
        %         %     legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southeast');
        %         %        f = gcf;
        %         %        saveas(f,sprintf('%s_spaces.fig',animal{kk}))
        %         %        exportgraphics(f,sprintf('Figures/%s_spaces.png',animal{kk}),'Resolution',300)
        %         %exportgraphics(f,sprintf('Figures/%s_spaces.fig',animal{kk}),'Resolution',300)
        %
        %         %        save(sprintf('ROI_%s.mat',animal{kk}),'ROI')
        %     end
    end
end

annotation('textbox',[.07 .88 .1 .1],'String','a','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.35 .88 .1 .1],'String','b','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.64 .88 .1 .1],'String','c','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
%annotation('textbox',[.5 .4 .1 .1],'String','d','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);

saveas(gcf,'fig_DR_S1.png')
saveas(gcf,'fig_DR_S1.fig')

f3 = figure('Position',  [100, 100, 800, 1150]);
%figure showing violin plots that replicate across DR strategy
for ii = 1:3
    for kk = 1:4
        
        endloop = 6;
        xlims = zeros(1,2);
        ylims = zeros(1,2);
        
        % kl is clusters
        for kl = 1:endloop
            dat = animals{kk}.clusters{ii}.data{kl};
            xlims = [xlims; [min(dat(:,1)) max(dat(:,1))]];
            ylims = [ylims; [min(dat(:,2)) max(dat(:,2))]];
        end
        xmax(kk) = max(abs(xlims(:)));
        ymax(kk) = max(abs(ylims(:)));
        
        % Make Grid over which to evaluate the RBF kernel density
        xlimits = [-1.1 1.1];
        ylimits = [-1.1 1.1];
        %xlimits = [floor(min(xlims(:,1)))-3 ceil(max(xlims(:,2)))+3];
        %ylimits = [floor(min(ylims(:,1)))-3 ceil(max(ylims(:,2)))+3];
        xpts = linspace(xlimits(1), xlimits(2),npts);
        ypts = linspace(ylimits(1), ylimits(2),npts);
        [X,Y] = meshgrid(xpts,ypts);
        xi = [X(:), Y(:)];
        XY{kk} = xi;
        
        for kl = 1:endloop
            
            hold on
            
            dat = animals{kk}.clusters{ii}.data{kl}./[xmax(kk) ymax(kk)];
            %plot_error_ellipse(animals{kk}.clusters{ii}.dat{kl},colors{kl},false)
            
            f = mvksdensity(dat,xi,'Bandwidth',kernel_width);
            f_reshape = reshape(f,size(X));
            % surf(X,Y,f_reshape)
            idx = find(f_reshape > prob_thresh*max(f));
            idxs{kl,kk} = idx;
            
            %h = plot(cloud{kl,ii}(rbf_boundary,1),cloud{kl,ii}(rbf_boundary,2),colors{kl},'LineWidth',2);
            %h = scatter(cloud(rbf_boundary,1),cloud(rbf_boundary,2),colors{kl},'LineWidth',2);
            %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            hold off
            ROI{kl,kk} = zeros(size(f_reshape));
            ROI{kl,kk}(idx) = 1;
            
        end
        whole_space = zeros(size(ROI{1,kk}));
        for kl = 1:endloop
            whole_space = whole_space | ROI{kl,kk};
            area(kl,kk) = sum(sum(ROI{kl,kk}));
        end
        whole_area(kk) = sum(sum(int64(whole_space)));
        
        % Plot whole space for sanity check
        %     figure()
        %     X = reshape(XY{kk}(:,1),npts,npts);
        %     Y = reshape(XY{kk}(:,2),npts,npts);
        %     bounds = bwboundaries(int64(whole_space));
        %     hold on
        %     for k = 1 : length(bounds)
        %         thisBoundary = bounds{k};
        %         x = thisBoundary(:, 2);
        %         y = thisBoundary(:, 1);
        %         vector = [];
        %         for l = 1:length(x)
        %             vector = [vector; [X(x(l),y(l)) Y(x(l),y(l))]];
        %         end
        %         h = plot(vector(:,2),vector(:,1),'Color',colors{kl},'Linewidth',2);
        %         if k >= 2
        %             h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %         end
        %     end
        %     xlabel('Dimension 1','FontSize',16)
        %     ylabel('Dimension 2','FontSize',16)
        %     title(sprintf('Whole Space, %s',animal{kk}),'FontSize',16)
        %     hold off
        
    end
    subplot(3,2,2*ii-1)
    normalized_area = area./whole_area;
    norm_area = normalized_area([6,2,1,5,3,4],:);
    labels = {spaces{6},spaces{2},spaces{1},spaces{5},spaces{3},spaces{4}};
    v = violinplot(norm_area',labels);
    ylim([0 1])
    title(sprintf('Normalized Area: %s',titles{ii}),'FontSize',16)
    ylabel('Fraction of total','FontSize',16)
    xtickangle(45)
    
    %starting the other side of the plots
    similarity = zeros(6,6,4);
    similarity_oneway = zeros(6,6,4);
    for kk = 1:4
        endloop = 6;
        for kl = 1:endloop
            for jl = 1:endloop
                %similarity matrix: param space x param space x animal
                % pulled from ROI so order is standard, sine... etc
                
                similarity(kl,jl,kk) = dice(ROI{kl,kk},ROI{jl,kk});
                %similarity_oneway(kl,jl,kk) = dice_one_way(ROI{kl,kk},ROI{jl,kk});
                %one-way dice: target is the columns
            end
        end
        separability(ii,kk) = norm(similarity(:,:,kk),'fro');
    end
    
    
    dice_avg = mean(similarity,3);
    dice_std = std(similarity,[],3);
    
    %loop to plot individual dice score graphs for each param space
    subplot(3,2,2*ii)
        
        tempinds = ones(1,6);
        tempinds(6) = 0;     %
        data_to_plot = transpose(squeeze(similarity(logical(tempinds),6,:)));
        
        violinplot(data_to_plot,labels);
        spaces_temp = spaces;
        set(gca,'XtickLabel',spaces_temp)
        xtickangle(45)
        ylabel('DICE Score','FontSize',16)
        ylim([0 1])
        title(sprintf('Pairwise similarity:  %s',titles{ii}),'FontSize',16)
    
end

annotation('textbox',[.05 .88 .1 .1],'String','a','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.05 .58 .1 .1],'String','b','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.05 .285 .1 .1],'String','c','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);

saveas(gcf,'fig_DR_S5.png')
saveas(gcf,'fig_DR_S5.fig')

f4 = figure('Position',  [100, 100, 800, 1150]);

% Makes all the Standard vs. Other space overlap plots
%kk = 1; % Use ARN088 spaces
first_space = [1, 4, 3, 2, 1]; % Standard, NPT, Standard
second_space = [6, 6, 6, 3, 2]; % Recording, Recording, Sine
for kk = 1:length(animals)
    for j = 1:length(first_space)
        subplot(length(first_space),4,length(animals)*j+kk-4)
        grp1 = first_space(j);
        grp2 = second_space(j);
        volume_ovl_2grps(j,kk) = sum(sum(ROI{grp1,kk} & ROI{grp2,kk}));
        dat1 = animals{kk}.clusters{ii}.data{grp1}./[xmax(kk) ymax(kk)];
        dat2 = animals{kk}.clusters{ii}.data{grp2}./[xmax(kk) ymax(kk)];
        intersect = double(ROI{grp1,kk}(:) & ROI{grp2,kk}(:));
        
        allpts = [dat1; dat2];
        
        vq = griddata(XY{kk}(:,1),XY{kk}(:,2),intersect,allpts(:,1),allpts(:,2),'nearest');
        pts = find(vq > 0);
        intersecting_pts = allpts(vq > 0,:);
        not_intersecting_pts = allpts(~(vq > 0),:);
        %scatter(not_intersecting_pts(:,1),not_intersecting_pts(:,2),sz,[0.6350, 0.0780, 0.1840],'filled')
        scatter(dat1(:,1),dat1(:,2),sz,colors{grp1},'filled');
        hold on
        scatter(dat2(:,1),dat2(:,2),sz,colors{grp2},'filled');
        %scatter(intersecting_pts(:,1),intersecting_pts(:,2),sz,[0, 0.5, 0],'filled')
        img = zeros(size(X));
        img(idxs{grp1,kk}) = 1;
        bounds = bwboundaries(img);
        for k = 1 : length(bounds)
            thisBoundary = bounds{k};
            x = thisBoundary(:, 2);
            y = thisBoundary(:, 1);
            vector = [];
            for l = 1:length(x)
                vector = [vector; [X(x(l),y(l)) Y(x(l),y(l))]];
            end
            h = plot(vector(:,2),vector(:,1),'Color',colors{grp1},'Linewidth',2,'HandleVisibility','off' );
            patch(vector(:,2), vector(:,1), colors{grp1},'FaceAlpha',0.4,'HandleVisibility','off');
            if k >= 2
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end
        
        img = zeros(size(X));
        img(idxs{grp2,kk}) = 1;
        bounds = bwboundaries(img);
        for k = 1 : length(bounds)
            thisBoundary = bounds{k};
            x = thisBoundary(:, 2);
            y = thisBoundary(:, 1);
            vector = [];
            for l = 1:length(x)
                vector = [vector; [X(x(l),y(l)) Y(x(l),y(l))]];
            end
            h = plot(vector(:,2),vector(:,1),'Color',colors{grp2},'Linewidth',2,'HandleVisibility','off');
            patch(vector(:,2), vector(:,1), colors{grp2},'FaceAlpha',0.4,'HandleVisibility','off');
            if k >= 2
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            
        end
        %         hold off
        %         xlabel('Dimension 1','FontSize',16)
        %         ylabel('Dimension 2','FontSize',16)
        if rem(length(animals)*j+kk-4,length(animals)) == 0
            legend({spaces{grp1},spaces{grp2}},'Location','Northwest')
        end
        
        if rem(length(animals)*j+kk-4,length(animals)) == 1
            ylabel({sprintf('%s vs.',spaces{grp1}),spaces{grp2}},'FontSize',14)
        end
        if (length(animals)*j+kk-4)<5
            title(sprintf('Subject %d',kk),'FontSize',14)
        end
        xlim(xlimits)
        ylim(ylimits)
        
        %f = gcf;
        %exportgraphics(f4,sprintf('Figures/%s_%s_vs_%s.png',animal{kk},spaces{grp1},spaces{grp2}),'Resolution',300)
    end
end
annotation('textbox',[.045 .88 .1 .1],'String','a','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.045 .7 .1 .1],'String','b','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.045 .525 .1 .1],'String','c','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.045 .345 .1 .1],'String','d','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.045 .185 .1 .1],'String','e','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);

saveas(f4,'fig_DR_S3.png')
saveas(f4,'fig_DR_S3.fig')


figure('Position',  [100, 100, 1150, 800])
subplot(1,2,1)
hold on
plot(mean_gamma,separability(1,:),'bo-','LineWidth',2);
plot(mean_gamma,separability(2,:),'ro-','LineWidth',2);
plot(mean_gamma,separability(3,:),'go-','LineWidth',2);
ylabel('Average DICE Score')
xlabel('Gamma Power Change (%)')
title('Max Amplitude Parameters')
legend(titles)
subplot(1,2,2)
hold on
plot(mean_gamma_low,separability(1,:),'bo-','LineWidth',2);
plot(mean_gamma_low,separability(2,:),'ro-','LineWidth',2);
plot(mean_gamma_low,separability(3,:),'go-','LineWidth',2);
ylabel('Average DICE Score')
xlabel('Gamma Power Change (%)')
legend(titles)
title('Min Amplitude Parameters')

saveas(gcf,'fig_DR_S4.png')
saveas(gcf,'fig_DR_S4.fig')
