%plot dimensionality-reduced PSD stuff for parameters idk

animal = {'ARN088','STV009','STV003','STV008'};
method = {'umap','tsne','nnmf','pca'};
subfield = {'reduction','coords','W','score'}; % MUST BE IN THIS ORDER, corresponds to methods
titles = {'UMAP','t-SNE','NNMF','PCA'};

scatter_type = 'log';
colors = {'g','m','c','r','b','k'};
sz = 10;

addpath('/Users/ERCOLE/Documents/Research/Repos/param_analysis/processed_data')

amp_thr_upper = [40,40,40,110];
amp_thr_lower = [20,20,20,70];
%amp_thr_upper = 110;
%amp_thr_lower = 70;

for kk = 2%1:length(animal)
    load(sprintf('%s_ParamAnalysis_DR.mat',animal{kk}));
    if strcmp(animal{kk},'STV003')
        endloop = 5;
    else
        endloop = 6;
    end
    figure()
    for ii = 1:length(method)
        subplot(2,2,ii)
        data = getfield(getfield(getfield(stim_dr_struct,method{ii}),scatter_type),subfield{ii});
        if ~strcmp(method{ii},'pca')
            scatter(data(:,1),data(:,2),sz,'filled');
        else
            scatter3(data(:,1),data(:,2),data(:,3),sz,'filled');
        end
        title(titles{ii})


        N_trials = stim_dr_struct.param_labels.N_trials;
        param_type = stim_dr_struct.param_labels.param_type;
    end
    sgtitle(sprintf('Raw Dimensionality Reduction, %s',animal{kk}))
    
    figure()
    for ii = 1:length(method)
        counter = 1;
        subplot(2,2,ii)
        data = getfield(getfield(getfield(stim_dr_struct,method{ii}),scatter_type),subfield{ii});
        if ~strcmp(method{ii},'pca')
            
            for kl = 1:endloop
                hold on
                scatter(data(counter:counter+N_trials(kl),1),data(counter:counter+N_trials(kl),2),sz,colors{kl},'filled');
                counter = counter + N_trials(kl);
                %disp(counter);
            end
            hold off
        else

            for kl = 1:endloop
                hold on
                scatter3(data(counter:counter+N_trials(kl),1),data(counter:counter+N_trials(kl),2),data(counter:counter+N_trials(kl),3),sz,colors{kl},'filled');
                counter = counter + N_trials(kl);
                %disp(counter);
            end
            hold off
        end
        title(titles{ii})
        legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southwest');
    end
    sgtitle(sprintf('Dimensionality Reduction, %s', animal{kk}))
    
    figure()
    for ii = 1:length(method)
        subplot(2,2,ii)
        data = getfield(getfield(getfield(stim_dr_struct,method{ii}),scatter_type),subfield{ii});
        counter = 1;
        for kl = 1:endloop
            if ~strcmp(stim_dr_struct.param_labels.param_type{kl},'Recording')
                param_temp = stim_dr_struct.param{kl};
                param_sel = param_temp(:,1)<amp_thr_upper(kk);
                points_to_plot = data(counter:counter+N_trials(kl),:);
                points_to_plot(param_sel,:) = [];
            else
                points_to_plot = data(counter:counter+N_trials(kl),:);
            end
            hold on
            if ~strcmp(method{ii},'pca')
                scatter(points_to_plot(:,1),points_to_plot(:,2),sz,colors{kl},'filled');
            else
                scatter3(points_to_plot(:,1),points_to_plot(:,2),points_to_plot(:,3),sz,colors{kl},'filled');
            end
            counter = counter + N_trials(kl);
            %disp(counter);
        end
        hold off
        legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southwest');
        title(titles{ii})
    end
    sgtitle(sprintf('High-Amplitude: %s',animal{kk}))
    
    figure()
    for ii = 1:length(method)
        subplot(2,2,ii)
        data = getfield(getfield(getfield(stim_dr_struct,method{ii}),scatter_type),subfield{ii});
        counter = 1;
        for kl = 1:endloop
            if ~strcmp(stim_dr_struct.param_labels.param_type{kl},'Recording')
                param_temp = stim_dr_struct.param{kl};
                param_sel = param_temp(:,1)>amp_thr_lower(kk);
                points_to_plot = data(counter:counter+N_trials(kl),:);
                points_to_plot(param_sel,:) = [];
            else
                points_to_plot = data(counter:counter+N_trials(kl),:);
            end
            hold on
            if ~strcmp(method{ii},'pca')
                scatter(points_to_plot(:,1),points_to_plot(:,2),sz,colors{kl},'filled');
            else
                scatter3(points_to_plot(:,1),points_to_plot(:,2),points_to_plot(:,3),sz,colors{kl},'filled');
            end
            counter = counter + N_trials(kl);
            %disp(counter);
        end
        hold off
        title(titles{ii})
        legend({'Standard','Sine','Sine (2N)','NPT','Poisson','Recording'},'Location','Southwest');
    end
    sgtitle(sprintf('Low-Amplitude: %s',animal{kk}))
end
