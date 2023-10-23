%param_dr_processing.m
%
%Processing script that concatenates PSD data for several different
%bandpower files corresponding to stimulation parameter spaces, performs
%several different dimensionality reduction strategies, and saves the
%results for each subject.

animal_IDs = {'STV009','STV003','ARN088','STV008'};

addpath(genpath('/Users/ERCOLE/Documents/Research/Repos/tools/umapFileExchange_2.2'));

for kk = 1:length(animal_IDs)
    tic
    cd('/Users/ERCOLE/Documents/Research/Data');
    
    load(sprintf('%s_Bandpower_Standard.mat',animal_IDs{kk}));
    bp_standard = bandpower_data;
    load(sprintf('%s_Bandpower_Sine.mat',animal_IDs{kk}));
    bp_sine = bandpower_data;
    load(sprintf('%s_Bandpower_Sine_2N.mat',animal_IDs{kk}));
    bp_sine2n = bandpower_data;
    load(sprintf('%s_Bandpower_NPT.mat',animal_IDs{kk}));
    bp_npt = bandpower_data;
    load(sprintf('%s_Bandpower_Poisson.mat',animal_IDs{kk}));
    bp_poisson = bandpower_data;
    
    count = 0;
    param_labels.trials = zeros(5,2);
    param_labels.param_type = cell(5,1);
    param_labels.N_trials = zeros(5,1);
    
    param_labels.N_trials(1) = bp_standard.config.N_trials;
    param_labels.param_type(1) = 'Standard';
    param_labels.trials(1,:) = [1 param_labels.N_trials(1)];
    param_labels.N_trials(2) = bp_sine.config.N_trials;
    param_labels.param_type(2) = 'Sine';
    param_labels.trials(2,:) = [param_labels.N_trials(1)+1 param_labels.N_trials(1)+param_labels.N_trials(2)];
    param_labels.N_trials(3) = bp_sine2n.config.N_trials;
    param_labels.param_type(3) = 'Sine_2N';
    param_labels.trials(3,:) = [param_labels.N_trials(2)+1 param_labels.N_trials(2)+param_labels.N_trials(3)];
    param_labels.N_trials(4) = bp_npt.config.N_trials;
    param_labels.param_type(4) = 'NPT';
    param_labels.trials(4,:) = [param_labels.N_trials(3)+1 param_labels.N_trials(3)+param_labels.N_trials(4)];
    param_labels.N_trials(5) = bp_poisson.config.N_trials;
    param_labels.param_type(5) = 'Poisson';
    param_labels.trials(5,:) = [param_labels.N_trials(4)+1 param_labels.N_trials(4)+param_labels.N_trials(5)];
    
    stim_psd = [bp_standard.stim_psd; bp_sine.stim_psd; bp_sine2n.stim_psd; bp_npt.stim_psd; bp_poisson.stim_psd];
    param = {bp_standard.param; bp_sine.param; bp_sine2n.param; bp_npt.param; bp_poisson.param};
    
    
    %unroll cell array into matrix for dimensionality reduction
    stim_psd = cell2mat(stim_psd);
        
    stim_psd_log = log10(stim_psd);
    
    [pca_coeff, pca_score, pca_latent, pca_tsquared, pca_explained, pca_mu] = pca(stim_psd);
    [pca_coeff_log, pca_score_log, pca_latent_log, pca_tsquared_log, pca_explained_log, pca_mu_log] = pca(stim_psd_log);
    
    [nnmf_w,nnmf_h,nnmf_d] = nnmf(stim_psd,2);
    [nnmf_w_log,nnmf_h_log_nnmf_d_log] = nnmf(stim_psd_log,2);
    
    rng(0);
    
    [tsne_coords,tsne_loss] = tsne(stim_psd);
    [tsne_coords_log,tsne_loss_log] = tsne(stim_psd_log);
    
    %NOW DO UMAP
    [umap_reduction, umap_class, umap_clusterIdentifiers, umap_extras]=run_umap(stim_psd);
    [umap_reduction_log, umap_class_log, umap_clusterIdentifiers_log, umap_extras_log]=run_umap(stim_psd_log);
    
    stim_dr_struct.param_labels = param_labels;
    stim_dr_struct.param = param;
    stim_dr_struct.stim_psd = stim_psd;
    
    stim_dr_struct.pca.coeff = pca_coeff;
    stim_dr_struct.pca.score = pca_score;
    stim_dr_struct.pca.latent = pca_latent;
    stim_dr_struct.pca.tsquared = pca_tsquared;
    stim_dr_struct.pca.explained = pca_explained;
    stim_dr_struct.pca.mu = pca_mu;
    
    stim_dr_struct.pca_log.coeff = pca_coeff_log;
    stim_dr_struct.pca_log.score = pca_score_log;
    stim_dr_struct.pca_log.latent = pca_latent_log;
    stim_dr_struct.pca_log.tsquared = pca_tsquared_log;
    stim_dr_struct.pca_log.explained = pca_explained_log;
    stim_dr_struct.pca_log.mu = pca_mu_log;
    
    stim_dr_struct.tsne.coords = tsne_coords;
    stim_dr_struct.tsne.loss = tsne_loss;
    stim_dr_struct.tsne_log.coords = tsne_coords_log;
    stim_dr_struct.tsne_log.loss = tsne_loss_log;
    
    stim_dr_struct.nnmf.W = nnmf_w;
    stim_dr_struct.nnmf.H = nnmf_h;
    stim_dr_struct.nnmf.D = nnmf_d;
    stim_dr_struct.nnmf_log.W = nnmf_w_log;
    stim_dr_struct.nnmf_log.H = nnmf_h_log;
    stim_dr_struct.nnmf_log.D = nnmf_d_log;
    
    %umap_reduction, umap_class, umap_clusterIdentifiers, umap_extras
    stim_dr_struct.umap.reduction = umap_reduction;
    stim_dr_struct.umap.class = umap_class;
    stim_dr_struct.umap.clusterIdentifiers = umap_clusterIdentifiers;
    stim_dr_struct.umap.extras = umap_extras; 
    stim_dr_struct.umap_log.reduction = umap_reduction_log; 
    stim_dr_struct.umap_log.class = umap_class_log; 
    stim_dr_struct.umap_log.clusterIdentifiers = umap_clusterIdentifiers_log;
    stim_dr_struct.umap_log.extras = umap_extras; 
    
    cd('/Users/ERCOLE/Documents/Research/Code/Repos/param_analysis/processed_data')
    save(sprintf('%s_ParamAnalysis_DR.mat',animal_IDs{kk}),stim_dr_struct);
    toc
    fprintf('%s file saved\n',animal_IDs{kk});
end