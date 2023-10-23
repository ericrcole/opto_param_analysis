%plot dimensionality-reduced PSD stuff for parameters idk

clear all;

% animal = {'ARN088','STV009','STV003','STV008'};
% method = {'umap','tsne','nnmf','pca'};
% subfield = {'reduction','coords','W','score'}; % MUST BE IN THIS ORDER, corresponds to methods
% titles = {'UMAP','t-SNE','NNMF','PCA'};

amp_sel = true;
amp_thr = [30 30 30 80];

animal = {'ARN088','STV009','STV003','STV008'};
method = {'umap','tsne','pca'};
subfield = {'reduction','coords','score'}; % MUST BE IN THIS ORDER, corresponds to methods
titles = {'UMAP','t-SNE','PCA'};

scatter_type = 'log';
%colors = {'g','m','c','r','b','k'};
sz = 10;

addpath('/Volumes/OPTO_PROC/OptoStim/dr_processed');
%addpath('~/Downloads/dr_processed');

amp_thr_upper = [40,40,40,110];
amp_thr_lower = [20,20,20,70];
%amp_thr_upper = 110;
%amp_thr_lower = 70;

for kk = 1:length(animal)
    load(sprintf('%s_ParamAnalysis_DR_3_nptadd.mat',animal{kk}));
    
    N_trials = stim_dr_struct.param_labels.N_trials;
    param_type = stim_dr_struct.param_labels.param_type;

    endloop = 6;

    for ii = 1:length(method)
        animals{kk}.animal = animal{kk};
        animals{kk}.clusters{ii}.volume = [];
        counter = 1;
        data = getfield(getfield(getfield(stim_dr_struct,method{ii}),scatter_type),subfield{ii});
        psd = stim_dr_struct.stim_psd;
        
        
        %if ~strcmp(method{ii},'pca')
        animals{kk}.clusters{ii}.method = method{ii};
        
        
%         xlims = zeros(1,2);
%         ylims = zeros(1,2);
% 
%         for kl = 1:endloop
%             xlims = [xlims; [min(data(:,1)) max(data(:,1))]];
%             ylims = [ylims; [min(data(:,2)) max(data(:,2))]];
%         end
%         xmax(ii) = max(abs(xlims(:)));
%         ymax(ii) = max(abs(ylims(:)));
% 
%         % Make Grid over which to evaluate the RBF kernel density
%         xlimits = [-1.1 1.1];
%         ylimits = [-1.1 1.1];
%         %xlimits = [floor(min(xlims(:,1)))-3 ceil(max(xlims(:,2)))+3];
%         %ylimits = [floor(min(ylims(:,1)))-3 ceil(max(ylims(:,2)))+3];
%         xpts = linspace(xlimits(1), xlimits(2),150);
%         ypts = linspace(ylimits(1), ylimits(2),150);
%         [X,Y] = meshgrid(xpts,ypts);
%         xi = [X(:), Y(:)];
%         XY{ii} = xi;
        startidx = 1;
        for kl = 1:endloop
            params = stim_dr_struct.param{kl};
            endidx = startidx+stim_dr_struct.param_labels.N_trials(kl)-1;
            dat = [data(startidx:endidx,1),data(startidx:endidx,2)];
            power_sd = psd(startidx:endidx,:);
            
            % Reject outliers whose minimum distance to other points in
            % the set is greater than the mean minimum distance + 3 std
            
            dist = sort(ipdm(dat));
            mindist = dist(2,:);
            m = mean(mindist);
            s = std(mindist);
            dat(mindist > m+3*s,:) = [];
            power_sd(mindist > m+3*s,:) = [];
            params(mindist > m+3*s,:) = [];
            
            if amp_sel && ~strcmp(param_type{kl},'Recording')
                param_exc = params(:,1)<amp_thr(kk);
                
                dat(param_exc,:) = [];
                power_sd(param_exc,:) = [];
                params(param_exc,:) = [];
            end
            
            animals{kk}.clusters{ii}.psd{kl} = power_sd;
            animals{kk}.clusters{ii}.param{kl} = params;
            
            if ~strcmp(method{ii},'nnmf')
%                 f = mvksdensity(dat,xi,'Bandwidth',0.05);
%                 f_reshape = reshape(f,size(X));
%                 % surf(X,Y,f_reshape)
%                 idx = find(f_reshape > 0.01*max(f));
%                 animals{kk}.clusters{ii}.idxs{kl} = idx;
%                 ROI = zeros(size(f_reshape));
%                 ROI(idx) = 1;
%                 animals{kk}.clusters{ii}.ROI{kl} = ROI;
                % Save cluster info
                
                animals{kk}.clusters{ii}.data{kl} = dat;
                [k,av] = convhull(dat(:,1),dat(:,2));
                animals{kk}.clusters{ii}.hull{kl} = k;
                animals{kk}.clusters{ii}.volume = [animals{kk}.clusters{ii}.volume av];
                animals{kk}.clusters{ii}.centroids{kl} = mean(dat);
                animals{kk}.clusters{ii}.covariance{kl} = cov(dat);
                [k,av] = boundary(dat);
                animals{kk}.clusters{ii}.boundary{kl} = k;
                animals{kk}.clusters{ii}.boundvol{kl} = av;
           
            end
%                 
            startidx = startidx + N_trials(kl);
            %disp(counter);
        end
        
    end

end

if amp_sel
    save('animals_v3_nptadd_30I.mat','animals')
else
    save('animals_v3_nptadd.mat','animals')
end