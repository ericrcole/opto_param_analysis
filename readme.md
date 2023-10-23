Code and processed data accompanying the paper "Irregular optogenetic stimulation waveforms can induce naturalistic patterns of hippocampal spectral activity"

Work in progress; will be finalized pending submission of the revised manuscript

Directory structure:

./exp_code: .m files for generating and plotting stimulation signals (as used in real-time TDT experiments)

./processing: code for processing raw TDT files into processed data files with segmented stimulation epochs, PSD and bandpower values, and generating neural latent space representations

./band_analysis: code for analysis of single-frequency biomarkers (Figs. 2 and 3)

./dr_analysis: code for analysis of dimension-reduced data (Figs. 4-6)

./data: processed data files.

Important files:

./dr_analysis/effect_of_stim_plots_v2.m: generates figure 4

./dr_analysis/plot_boundary_size_fig_v3.m: generates figure 5

./dr_analysis/plot_query_pts_v5_erc.m: generates subplots of figure 6


Data structure:

Each file represents the data collected for one subject and for one parameter space, concatenated for all separate experiment sessions. Each file is a .mat data structure with several fields:
- fields named after canonical bands contain an Nx2 matrix, where N is the number of stimulation trials, column 1 is the pre-stimulation value, and column 2 is the during-stimulation value
- each row of .param contains the stimulation parameters applied at the given trial
