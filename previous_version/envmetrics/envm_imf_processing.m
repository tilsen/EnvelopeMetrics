
%envelope parameters
par.bandpass = [400 4000];  %: vocalic energy bandpass range (there is no a priori justification for this range)
par.lowpass = 10;           %: lowpass filtering cutoff
par.ds = 100;               %: downsampling factor

%% initial extraction of instantaneous frequencies:
files = {'example.wav','example2.wav'};
for i=1:length(files)
    [x,par.Fs] = audioread(files{i}); 
    env = envm_band_energy(x,par);
    
    %normalize:
    env = env-mean(env);
    env = env/max(abs(env));  
    
    %window
    envw = tukeywin(length(env),0.2).*env;
    
    %empirical mode decomposition (matlab built-in)
    imf = emd(envw,'SiftRelativeTolerance',0.1,'MaxNumIMF',4);
    
    %hilbert-huang transform
    envelope_Fs = par.Fs/par.ds;
    [~,~,~,wi] = hht(imf,envelope_Fs);
    
    IMF{i} = num2cell(imf,1);
    W{i} = num2cell(wi,1);
end

%1. the first and last 100ms (approximately) can exhibit edge artifacts, so
%replace them with NaN:
edgedur = 0.1;

for i=1:length(W) %loop over chunks
    for j=1:length(W{i}) %loop over IMFs
        W{i}{j}(1:round(envelope_Fs*edgedur)) = nan;
        W{i}{j}(end-round(envelope_Fs*edgedur)+1:end) = nan;
    end
end

%2. impose bounds on frequencies
lower_bound = 0;
upper_bound = 13.16; %frequency at which magnitude response of 4th order butterworth is -10dB
for i=1:length(W) %loop over chunks
    for j=1:length(W{i}) %loop over IMFs
        W{i}{j}(W{i}{j}<lower_bound) = nan;
        W{i}{j}(W{i}{j}>upper_bound) = nan;
    end
end

%3. for each IMF, exclude values above some percentile

exc_perc = 99;

%first calculate the exclusion thresholds
for i=1:2 %loop over IMF #s
    allw = cellfun(@(c)c{i}',W,'un',0);
    allw = [allw{:}];
    wthresh(i) = prctile(allw,exc_perc);
end

%now replace with nan with percentiles
for i=1:length(W) %loop over chunks
    for j=1:2 %loop over IMFs
        W{i}{j}(W{i}{j}>wthresh(j)) = nan;
    end
end




