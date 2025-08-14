function [pars] = envm_default_params()

defPars.bandpass = [400 nan]; %upper band edge is determined from Fs
defPars.lowpass = 10;
defPars.ds = 100; %downsampling factor
defPars.tukeywin_param = 0.1;
defPars.nfft = 2048; 
defPars.sm_Hz = 1;
defPars.powerratio_freq_bins = [1 3.5 10];
defPars.centroid_freq_bins = [1 10];
defPars.max_imf = 3;
defPars.edge_null = 0.1;
defPars.imf_freq_bounds = [0 13.16]; %frequency at which magnitude response of 4th order butterworth is -10dB
defPars.sift_relative_tol = 0.1;
defPars.imf_freq_exclusion_prctile = 99;
defPars.verbose = false;

descriptions = {
    'bandpass filter edges. NaN = determine from sampling rate.'
    'lowpass filter cutoff'
    'downsampling factor' 
    'tukey window parameter for envelopes'
    'number of dft coefficients'
    'power spectrum smoothing bandwidth'
    'frequency bins for LFFA power ratios'
    'frequency range for LFFA centroid'
    'maximum number of IMFs to return'
    'period of imf margins to ignore'
    'allowed range of imf frequencies'
    'sift relative tolerance for Matlab built-in emd function'
    'percentile for exclusion of imf frequencies across dataset'
    'verbose output'
    };

ff = fieldnames(defPars);

pars = table(ff,'VariableNames',{'paramName'});

for i=1:height(pars)
    pars.value{i} = defPars.(pars.paramName{i});
end

pars.description = descriptions;


end