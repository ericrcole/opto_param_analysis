
load('animals.mat')
addpath('Violinplot-Matlab')
animal = {'ARN088','STV009','STV003','STV008'};
ii = 1; % Use UMAP
plot_animal = 1; % only plot ARN088
titles = {'UMAP','t-SNE','PCA'};

spaces = {'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'};

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

for kk = 1:4
    if strcmp(animal{kk},'STV003')
        endloop = 5;
    else
        endloop = 6;
    end
    
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
        figure()
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
        title('Dimensionality Reduced Space','FontSize',16)
        xlabel('Dimension 1','FontSize',16)
        ylabel('Dimension 2','FontSize',16)
        legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southeast');
        xlim(xlimits)
        ylim(ylimits)
        f = gcf;
        exportgraphics(f,sprintf('Figures/%s_data.png',animal{kk}),'Resolution',300)
    end

    if kk == plot_animal
        figure()
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
        legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southeast');
        xlabel('Dimension 1','FontSize',16)
        ylabel('Dimension 2','FontSize',16)
        title(sprintf('Boundaries of Parameter Spaces'),'FontSize',16)
        xlim(xlimits)
        ylim(ylimits)
        hold off
%     legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southeast');
        f = gcf;
        exportgraphics(f,sprintf('Figures/%s_spaces.png',animal{kk}),'Resolution',300)
        save(sprintf('ROI_%s.mat',animal{kk}),'ROI')
    end
end

%% Analyze ROIs

% Calculate total DR space volume in pixels for each animal
for kk = 1:4
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

figure()
normalized_area = area./whole_area;
norm_area = normalized_area([6,2,1,5,3,4],:);
labels = {spaces{6},spaces{2},spaces{1},spaces{5},spaces{3},spaces{4}};
v = violinplot(norm_area',labels);
title('Normalized area of stimulation spaces','FontSize',16)
ylabel('Percent of total area','FontSize',16)
f = gcf;
exportgraphics(f,'Figures/normalized_area.png','Resolution',300)

% Makes all the Standard vs. Other space overlap plots
kk = 1; % Use ARN088 spaces
first_space = [1, 4, 1]; % Standard, NPT, Standard
second_space = [6, 6, 2]; % Recording, Recording, Sine
for j = 1:3
    figure()
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
    scatter(not_intersecting_pts(:,1),not_intersecting_pts(:,2),sz,[0.6350, 0.0780, 0.1840],'filled')
    hold on
    scatter(intersecting_pts(:,1),intersecting_pts(:,2),sz,[0, 0.5, 0],'filled')
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
        h = plot(vector(:,2),vector(:,1),'Color',colors{grp1},'Linewidth',2);
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
        h = plot(vector(:,2),vector(:,1),'Color',colors{grp2},'Linewidth',2);
        if k >= 2
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
    end
    hold off
    xlabel('Dimension 1','FontSize',16)
    ylabel('Dimension 2','FontSize',16)
    legend({'Non-overlapping','Overlapping',spaces{grp1},spaces{grp2}},'Location','Northeast')
    title(sprintf('%s vs. %s',spaces{grp1},spaces{grp2}),'FontSize',16)
    xlim(xlimits)
    ylim(ylimits)
    f = gcf;
    exportgraphics(f,sprintf('Figures/%s_%s_vs_%s.png',animal{kk},spaces{grp1},spaces{grp2}),'Resolution',300)
end


% Calculate Pair-wise set similarity with Dice Score
similarity = zeros(6,6,4);
for kk = 1:4
    for kl = 1:endloop
        for jl = 1:endloop
            similarity(kl,jl,kk) = dice(ROI{kl,kk},ROI{jl,kk});
        end
    end
end
    
dice_avg = mean(similarity,3);
dice_std = std(similarity,[],3);




