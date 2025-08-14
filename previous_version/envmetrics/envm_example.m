
%% A set of audio chunks can be batch-processed as follows:
files = dir('*.wav');
X={};
for i=1:length(files)
	[X{i},Fs] = audioread(files(i).name);
end
par.Fs = Fs;

metrics = envm_metrics_batch(X,par,'verbose',true);
disp(metrics);


%% Here is how to extract metrics from a single file:

[wav,Fs] = audioread('example.wav');
t_wav = linspace(0,length(wav)-1,length(wav))/Fs;

par.bandpass = [400 4000];  %: vocalic energy bandpass range (there is no a priori justification for this range)
par.lowpass = 10;           %: lowpass filtering cutoff
par.Fs = Fs;                %: sampling rate
par.ds = 100;               %: downsampling factor

%extract the envelope:
[env,t_env] = envm_band_energy(wav,par);

%normalize:
env = env-mean(env);
env = env/max(abs(env));

%tukey windowing:
envw = tukeywin(length(env),0.2).*env;

%smoothed power spectrum
par.nfft = 2048; 
par.sm_Hz = 1;
par.Fs = par.Fs/par.ds;

[smpsd,f] = envm_smoothed_psd(envw,par);

%% emd/hht
%empirical mode decomposition/hht (using a modified version of 
%Alan Tan's toolbox (http://www.mathworks.com/matlabcentral/fileexchange/19681-hilbert-huang-transform/content/plot_hht.m) 

%[imf,w,t_imf] = envm_hht(envw,1/par.Fs); %takes sampling period instead of freq.

%using matlab built-in functions (close but not quite identical to the toolbox):
imf = emd(envw,'SiftRelativeTolerance',0.1,'MaxNumIMF',4);
[hs,~,t_imf,w] = hht(imf,par.Fs);
imf = num2cell(imf,1);
w = num2cell(w,1);

%% post-processing IMF instantaneous frequencies
%the IMF frequencies diverge so it is important to exclude ones that are
%outside of meaningful ranges. This should be done using an entire corpus of chunks. see envm_imf_processing for examples


%% plot waveform, envelope, smoothed envelope spectrum, imfs, and instantaneous frequencies
figure;
subplot(4,3,1); plot(t_wav,wav); axis tight;
subplot(4,3,4); plot(t_env,envw); axis tight;
subplot(4,3,[7 10]); plot(f,smpsd); 
hold on; axis tight; title('smoothed power spectrum'); xlabel('frequency');

for i=1:length(imf)
    subplot(4,3,3*i-1); plot(t_imf,imf{i}); hold on; axis tight; title(sprintf('imf %i',i));
    subplot(4,3,3*i); plot(t_imf,w{i}); hold on; axis tight; title(sprintf('imf %i instantaneous freq.',i));
end

%% metrics
par.powerratio_freq_bins = [1 3.5 10];
par.centroid_freq_bins = [1 10];

SPEC = envm_psd_metrics(smpsd,f,par) %#ok<*NOPTS>
EMD = envm_emd_metrics(imf,w)




