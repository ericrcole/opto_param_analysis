%produces final version of figure 4 and one of the supplements

save_subfigs = false;
version_string = '_v3_nptadd_30I';
load(sprintf('animals%s.mat',version_string))
addpath('Violinplot-Matlab')
animal = {'ARN088','STV009','STV003','STV008'};
ii = 1; % Use UMAP
plot_animal = 1; % only plot ARN088
titles = {'UMAP','t-SNE','PCA'};

spaces = {'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'};

% col_scale = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660, 0.6740, 0.1880],'black'};
% labels = {'Sine','Standard','Poisson','Nested Pulse','Double Sine','Behavior'}

scatter_type = 'log';
%colors = {'g','m','c','r','b','k'};
colors = {[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.4660, 0.6740, 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],'k'};
%colors = {[0 0.30 0.7810],[0.8500 0.5 0.0980],[0.4940 0.1840 0.5560],[0.4660 0.9 0.1880],[0.9 0.0780 0.1840],'k'};
nic = [0.6350, 0.0780, 0.1840];
ic = [0, 0.5, 0];
sz = 10;

npts = 200;
kernel_width = 0.05; %0.05
prob_thresh = 0.005;  %0.01


xmax = zeros(3,1);

f1 = figure('Position',  [100, 100, 1150, 800]);

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
    if kk == plot_animal
        %figure()
    end
    for kl = 1:endloop
        
        hold on
       
        dat = animals{kk}.clusters{ii}.data{kl}./[xmax(kk) ymax(kk)];
        %plot_error_ellipse(animals{kk}.clusters{ii}.dat{kl},colors{kl},false)
         
        f = mvksdensity(dat,xi,'Bandwidth',kernel_width);
        f_reshape = reshape(f,size(X));
        % surf(X,Y,f_reshape)
        idx = find(f_reshape > prob_thresh*max(f));
        idxs{kl,kk} = idx;
        
        if kk == plot_animal
             subplot(2,3,1)
            scatter(dat(:,1),dat(:,2),sz,colors{kl},'filled')
        end
        
        %h = plot(cloud{kl,ii}(rbf_boundary,1),cloud{kl,ii}(rbf_boundary,2),colors{kl},'LineWidth',2);
        %h = scatter(cloud(rbf_boundary,1),cloud(rbf_boundary,2),colors{kl},'LineWidth',2);
        %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        hold off
        ROI{kl,kk} = zeros(size(f_reshape));
        ROI{kl,kk}(idx) = 1;
        
    end
    if kk == plot_animal
       
        title('Neural latent space','FontSize',16)
        xlabel('Dimension 1','FontSize',16)
        ylabel('Dimension 2','FontSize',16)
        legend({'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'},'Location','Southeast');
        xlim(xlimits)
        ylim(ylimits)
        f = gcf;
        if save_subfigs
        saveas(f,sprintf('%s_DR.fig',animal{kk}));
        exportgraphics(f,sprintf('Figures/%s_data.png',animal{kk}),'Resolution',300)
        %exportgraphics(f,sprintf('Figures/%s_data.fig',animal{kk}),'Resolution',300)
        end
    end

    if kk == plot_animal
        subplot(2,3,2)
        for kl = 1:endloop
            hold on
            dat = animals{kk}.clusters{ii}.data{kl}./[xmax(kk) ymax(kk)];
            img = zeros(size(X));
            img(idxs{kl,kk}) = 1;
            bounds = bwboundaries(img);
            
            for k = 1 : length(bounds)
                thisBoundary = bounds{k};
                x = thisBoundary(:, 2);
                y = thisBoundary(:, 1);
                vector = [];
                for l = 1:length(x)
                    vector = [vector; [X(x(l),y(l)) Y(x(l),y(l))]];
                end
                h = plot(vector(:,2),vector(:,1),'Color',colors{kl},'Linewidth',2);
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            scatter(dat(:,1),dat(:,2),sz,colors{kl},'filled')
        end
        legend({'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'},'Location','Southeast');
        xlabel('Dimension 1','FontSize',16)
        ylabel('Dimension 2','FontSize',16)
        title(sprintf('Parameter space boundaries'),'FontSize',16)
        xlim(xlimits)
        ylim(ylimits)
        hold off
%     legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southeast');
%        f = gcf;
%        saveas(f,sprintf('%s_spaces.fig',animal{kk}))
%        exportgraphics(f,sprintf('Figures/%s_spaces.png',animal{kk}),'Resolution',300)
        %exportgraphics(f,sprintf('Figures/%s_spaces.fig',animal{kk}),'Resolution',300)
        
%        save(sprintf('ROI_%s.mat',animal{kk}),'ROI')
    end
end

%% Analyze ROIs

% Calculate total DR space volume in pixels for each animal
for kk = 1:4
    
        endloop = 6;
    
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

%%
subplot(2,3,3)
normalized_area = area./whole_area;
norm_area = normalized_area([6,2,1,5,3,4],:)';
%labels = {spaces{6},spaces{2},spaces{1},spaces{5},spaces{3},spaces{4}};
%labels = {'Standard','Sine','Double Sine','Nested Pulse','Poisson','Behavior'};
labels = {'Behavior','Sine (2D)','Standard (3D)','Poisson (3D)','Double Sine (4D)', 'Nested Pulse (5D)'};
%v = violinplot(norm_area',labels);
plot_ord = [6,2,1,5,3,4];
for kk = 1:length(plot_ord)
    if kk == 1
        Violin(norm_area(:,kk),kk,'ViolinColor',[0 0 0]);
        
        text(kk,0.9,'***','FontSize',28,'HorizontalAlignment','center')
    else
        Violin(norm_area(:,kk),kk,'ViolinColor',colors{plot_ord(kk)});
    end
end
ylim([0 1])
xticks(1:length(plot_ord))
xticklabels(labels)
title('Normalized area','FontSize',16)
ylabel('Fraction of total','FontSize',16)
xtickangle(45)
set(gca,'fontsize',16)

f = gcf;
%exportgraphics(f,'Figures/normalized_area.png','Resolution',300)

% Makes all the Standard vs. Other space overlap plots
kk = 1; % Use ARN088 spaces
first_space = [1, 4, 1]; % Standard, NPT, Standard
second_space = [6, 6, 2]; % Recording, Recording, Sine
for j = 1:2
    subplot(2,3,3+j)
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
    rec_bounds = bounds;
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
    hold off
    xlabel('Dimension 1','FontSize',16)
    ylabel('Dimension 2','FontSize',16)
    legend({spaces{grp1},spaces{grp2}},'Location','Southeast')
    title(sprintf('%s vs. %s',spaces{grp1},spaces{grp2}),'FontSize',16)
    xlim(xlimits)
    ylim(ylimits)
    f = gcf;
    %exportgraphics(f,sprintf('Figures/%s_%s_vs_%s.png',animal{kk},spaces{grp1},spaces{grp2}),'Resolution',300)
end

%%
%saveas(gcf,'param_boundary_size.png')
%saveas(gcf,'param_boundary_size.fig')

% Calculate Pair-wise set similarity with Dice Score
similarity = zeros(6,6,4);
similarity_oneway = zeros(6,6,4);
for kk = 1:4
    endloop = 6;
    for kl = 1:endloop
        for jl = 1:endloop
            %similarity matrix: param space x param space x animal
                % pulled from ROI so order is standard, sine... etc
            
            similarity(kl,jl,kk) = dice(ROI{kl,kk},ROI{jl,kk});
            similarity_oneway(kl,jl,kk) = dice_one_way(ROI{kl,kk},ROI{jl,kk});
            %one-way dice: target is the columns
        end
    end
end

dice_avg = mean(similarity,3);
dice_std = std(similarity,[],3);

dice_ow_avg = mean(similarity_oneway,3);
dice_ow_std = std(similarity_oneway,[],3);

%loop to plot individual dice score graphs for each param space
subplot(2,3,6)
for k1 = 6:6
    subplot(2,3,k1)
%     tempinds = ones(1,6); 
%     tempinds(k1)=0;
    data_to_plot = transpose(squeeze(similarity(:,k1,:)));
    data_to_plot(:,k1) = [];
    
    %violinplot(data_to_plot,labels);
    for kk = 1:size(data_to_plot,2)
        if kk == 1
            text(kk,0.9,'***','FontSize',28,'HorizontalAlignment','center')
        end
        Violin(data_to_plot(:,kk),kk,'ViolinColor',colors{kk});
    end
    ylim([0 1])
    xticks(1:length(plot_ord))
    xticklabels(labels)

    spaces_temp = spaces;
    spaces_temp(k1) = [];
    
    set(gca,'XtickLabel',spaces_temp)
    set(gca,'fontsize',16)
    xtickangle(45)
    ylabel('DICE Score','FontSize',16)
    title('Pairwise similarity vs. Behavior','FontSize',16)
end

annotation('textbox',[.08 .88 .1 .1],'String','a','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.36 .88 .1 .1],'String','b','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.645 .88 .1 .1],'String','c','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.08 .4 .1 .1],'String','d','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.36 .4 .1 .1],'String','e','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.645 .4 .1 .1],'String','f','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);



saveas(gcf,'fig4-2.png')
saveas(gcf,'fig4-2.fig')

f2 = figure('Position',  [100, 100, 850, 600]);
for k1 = 1:6
    subplot(2,3,k1)
    tempinds = ones(1,6); 
    tempinds(k1)=0;
    data_to_plot = transpose(squeeze(similarity(logical(tempinds),k1,:)));
        
    violinplot(data_to_plot,labels);
    spaces_temp = spaces;
    spaces_temp(k1) = [];
    set(gca,'XtickLabel',spaces_temp,'FontSize',14)
    xtickangle(45)
    ylabel('DICE Score','FontSize',14)
    ylim([0 1])
    title(spaces{k1},'FontSize',16)
    
    
end

annotation('textbox',[.08 .88 .1 .1],'String','a','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.36 .88 .1 .1],'String','b','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.645 .88 .1 .1],'String','c','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.08 .4 .1 .1],'String','d','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.36 .4 .1 .1],'String','e','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);
annotation('textbox',[.645 .4 .1 .1],'String','f','FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 24);


saveas(gcf,'param_dice2.png')
saveas(gcf,'param_dice2.fig')

f3 = figure('Position',  [100, 100, 1150, 800]);
for k1 = 1:6
    subplot(2,3,k1)
    tempinds = ones(1,6); 
    tempinds(k1)=0;
    data_to_plot = transpose(squeeze(similarity_oneway(logical(tempinds),k1,:)));
        %this is right, just trust me - each subplot is the param space
        %being projected *onto*
       
    violinplot(data_to_plot,labels);
    spaces_temp = spaces;
    spaces_temp(k1) = [];
    set(gca,'XtickLabel',spaces_temp)
    xtickangle(45)
    ylabel('Adjusted DICE')
    title(spaces{k1})
end
suptitle('Param Space Projection Coverage')

saveas(gcf,'param_dice_one_way2.png')
saveas(gcf,'param_dice_one_way2.fig')

%% block to run/test stats
[p1,tbl,stats]=friedman(norm_area',1);
mcomp = multcompare(stats);
disp(p1)

[p2,tbl,stats]=friedman(data_to_plot,1);
mcomp = multcompare(stats);
disp(p2)

%% bootstrapped statistics
%sample data with replacement and re-compute stats



% area_std = norm_area(2,:);
% area_other = norm_area(3:end,:); area_other = area_other(:);
% [p,h,stats] = ranksum(area_std,area_other);
% disp(p)
% disp(h)
% disp(stats)
% 
% temp = padcat(area_std',area_other);
% [p,tbl,stats] = anova1(temp);
% 
% [p,tbl,stats] = anova1(norm_area');
% % tbl = array2table(results,"VariableNames", ...
% %     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% res = multcompare(stats);
% res = array2table(res,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% 
% data_to_plot = transpose(squeeze(similarity(:,k1,:)));
% data_to_plot(:,k1) = [];
% 
% [p,tbl,stats] = anova1(data_to_plot);
% % tbl = array2table(results,"VariableNames", ...
% %     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% res = multcompare(stats);
% res = array2table(res,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% 
% temp = data_to_plot(2,:); temp = temp(:);
% temp = padcat(data_to_plot(:,1),temp);
% [p,tbl,stats] = anova1(temp);
% 
% %ranksum test
% temp = data_to_plot(:,2:end); temp = temp(:);
% temp = padcat(data_to_plot(:,1),temp);
% [p,h,ci] = ranksum(temp(:,1),temp(:,2));
% %signrank test
% tempmat = data_to_plot(:,2:end);
% diffs_sr = tempmat - data_to_plot(:,1);
% [p,h,ci] = signrank(diffs_sr(:));
% 
% %repeated measures anova
% % rm_arr1 = num2cell(data_to_plot(:));
% % rm_arr2 = [repmat(labels(2),4,1); repmat(labels(3),4,1);repmat(labels(4),4,1);repmat(labels(5),4,1);repmat(labels(6),4,1)];
% % rm_arr3 = repmat({'A';'B';'C';'D'},5,1);
% % rm_arr = [rm_arr1,rm_arr2,rm_arr3];
% % rm_arr = cell2table(rm_arr,"VariableNames",["Value","Param","Subject"]);
% 
% %rm_arr = array2table(data_to_plot',"VariableNames",{'S1','S2','S3','S4'});
% %rm_arr.param = labels([2 1 3 4 5])';
% % Meas = table([1 2 3 4]','VariableNames',"Subject");
% % rm = fitrm(rm_arr,'S1-S4 ~ param','WithinDesign',Meas);
% 
% rm_arr = array2table(data_to_plot,"VariableNames",{'Std','Sine','Psn','DblSne','NPT'});
% rm_arr.subject = {'A','B','C','D'}';
% 
% Meas = table(categorical([1 2 3 4]'),'VariableNames',"Subject");
% rm = fitrm(rm_arr,'Std-NPT ~ 1','WithinDesign',Meas);


function dice = dice_one_way(roi1,roi2)
    % one-way dice = intersection of roi1, roi2 divided by area of roi2
    dice = sum(roi1 & roi2)/sum(roi2);
end

function plot_patch(x,mean_tr,dev_tr,col)
    
x_new = [x fliplr(x)];
y_plot = [mean_tr - dev_tr, fliplr(mean_tr + dev_tr)];
patch(x_new, y_plot, col,'FaceAlpha',0.7);

end