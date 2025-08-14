function [be,t,bp] = envm_band_energy(wav,par)

%INPUTS
% wav: the acoustic signal
% par.bandpass: vocalic energy bandpass range 
% par.lowpass: lowpass filtering cutoff
% par.Fs: sampling rate
% par.ds: downsampling factor

%OUTPUTS
%be: band energy (e.g. vocalic energy  amplitude envelope), unnormalized
%t:  time vector
%bp: band-pass filtered signal

if ~isfield(par,'bandpass_order'), par.bandpass_order = 4; end
if ~isfield(par,'lowpass_order'), par.lowpass_order = 4; end

%PARAMETERS
Fs              = par.Fs;
ds              = par.ds;     
bandpass        = par.bandpass;                           %bandpass frequencies
lowpass         = par.lowpass;                            %lowpass cutoff frequency
nyquist         = Fs/2;                                     %nyquist frequency
bandpass_norm   = bandpass/nyquist;                         %normalized bandpass vector
bandpass_order  = par.bandpass_order;
lowpass_norm    = lowpass/nyquist;                          %normalized lowpass cutoff
lowpass_order   = par.lowpass_order;

%PROCESSING
wav_dc          = wav - mean(wav);                          %remove DC
[bbp,abp]       = butter(bandpass_order,bandpass_norm);     %1st order, bandpass filter pareters
sig_bp          = filtfilt(bbp,abp,wav_dc);                   %bandpass filter
%sig_power      = (abs(sig_bp).^2);                         %power
sig_mag         = (abs(sig_bp));                            %magnitude
[blp,alp]       = butter(lowpass_order,lowpass_norm);                   %4th order, lowpass (<20Hz) filter pareters
sig_lp          = filtfilt(blp,alp,sig_mag);                  %lowpass filter
sig_ds          = downsample(sig_lp,ds);                   %downsample
sig_dc          = (sig_ds - mean(sig_ds));                  %remove DC


%calculate phase_delay
% [phi,w,s]=phasedelay(bbp,abp,N,Fs);
% inds = find(w>bandpass(1) & w<bandpass(2));
% delay_bp = mean(phi(inds))/Fs;
% [phi,w,s]=phasedelay(blp,alp,N,Fs);
% inds = find(w<lowpass);
% delay_lp = mean(phi(w>0 & w<lowpass))/Fs;
% phase_delay = delay_bp + delay_lp;


%OUTPUTS
be = sig_dc;
bp = sig_bp;
t = linspace(0,length(be)-1,length(be))/(Fs/ds);

end