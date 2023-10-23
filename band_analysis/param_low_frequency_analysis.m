%param_low_frequency_analysis.m
%
%script to analyze differences in driving low-frequency activity as a
%result of multiple different stimulation parameters and parameter spaces.

animal = {'ARN088','STV009','STV008','STV003'};

cd('/Users/ERCOLE/Documents/Research/Data');

for kk = 1:length(animal)
    
load(sprintf('%s_Bandpower_Awake_Standard.mat',animal{kk}));
bp_standard = bandpower_data;
load(sprintf('%s_Bandpower_Awake_Sine.mat',animal{kk}));
bp_sine = bandpower_data;
load(sprintf('%s_Bandpower_Awake_Sine_2N.mat',animal{kk}));
bp_sine2n = bandpower_data;
load(sprintf('%s_Bandpower_Awake_NPT.mat',animal{kk}));
bp_npt = bandpower_data;
load(sprintf('%s_Bandpower_Awake_Poisson.mat',animal{kk}));
bp_poisson = bandpower_data;

figure
%difference in driving oscillations between standard pulse and
    %sinusoidal stimulation
    
    
end


    %difference in driving oscillations between standard pulse and
    %sinusoidal stimulation
