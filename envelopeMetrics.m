classdef envelopeMetrics
    
    properties
        X
        Fs      
        env_Fs
        env_Passband = [400 4000]
        env_Lowpass = 10
        env_BandpassFilterOrder = 4
        env_LowpassFilterOrder = 4
        env_Downsample = 100
        env_Rescale = true        
        env_TukeywinParam = nan
        env_EdgeAttenutation = 0.05 %period of time to attenuate edges
        spec_Nfft = 2048
        spec_SmoothBw = 1
        spec_PowerBins = [1 3.5; 3.5 10]
        spec_CentroidBins = [1 10]        
        emd_MaxImf = 3
        emd_EdgeNull = 0.1
        emd_ImfFreqBounds = [0 13.16] %frequency at which magnitude response of 4th order butterworth is -10dB
        emd_SiftRelTol = 0.1
        emd_FreqExclusionPercentile = 99
        verbose = false
    end

    methods
        function obj = envelopeMetrics(X,Fs)
            if ~iscell(X) || ~(any(size(X{1})==1))
                fprintf('ERROR: input waveforms as a cell array of row or column vectors');
                return
            end
            if nargin==1 || ~isscalar(Fs)
                fprintf('ERROR: sampling rate must be provided as second input');
                return
            end            
            obj.X = cellfun(@(c){c(:)'},X);
            obj.Fs = Fs;
            obj.env_Fs = obj.Fs/obj.env_Downsample;
        end

        function [metrics,envelopes,times,spectra,freqs,imfs,imfw] = getMetrics(obj)
            [envelopes,times,~,~,obj] = obj.extractEnvelopes();
            [spectra, freqs] = obj.extractSpectra(envelopes);
            spectralMetrics = obj.spectralMetrics(spectra,freqs);
            [emdMetrics,imfs,imfw] = obj.emdMetrics(envelopes);
            metrics = [spectralMetrics emdMetrics];
            cols = metrics.Properties.VariableNames;
            colsOrdered = cols(~cellfun('isempty',(regexp(cols,'sbpr|scntr|imf_ratio', 'once'))));
            colsOrdered = [colsOrdered sort(setdiff(cols,colsOrdered))];
            metrics = metrics(:,colsOrdered);
        end

        function [envelopes, times,passbandFiltered, lowpassFiltered,obj] = extractEnvelopes(obj)            
            x = cellfun(@(c){c-mean(c)},obj.X);

            if isnan(obj.env_Passband)
                obj.env_Passband(2) = obj.Fs/2;
            end

            [bbp,abp] = butter(obj.env_BandpassFilterOrder,obj.env_Passband/(obj.Fs/2));
            [blp,alp] = butter(obj.env_LowpassFilterOrder,obj.env_Lowpass/(obj.Fs/2));
            
            passbandFiltered = cellfun(@(c){filtfilt(bbp,abp,c)},x);
            lowpassFiltered = cellfun(@(c){filtfilt(blp,alp,abs(c))},passbandFiltered);
            
            envelopes = cellfun(@(c){downsample(c,obj.env_Downsample)},lowpassFiltered);
            obj.env_Fs = obj.Fs/obj.env_Downsample;

            if obj.env_Rescale          
                envelopes = cellfun(@(c){c/max(abs(c))},envelopes);
            end
            times = cellfun(@(c){(0:length(c)-1)/obj.env_Fs},envelopes);

        end

        function [envelopes] = attenuateEdges(obj,envelopes)
            if ~isempty(obj.env_TukeywinParam) && ~isnan(obj.env_TukeywinParam)
                envelopes = cellfun(@(c){c.*tukeywin(length(c),obj.env_TukeywinParam)'},envelopes);
            end
            if ~isempty(obj.env_EdgeAttenutation) && ~isnan(obj.env_EdgeAttenutation)
                ec = obj.env_EdgeAttenutation;
                times = cellfun(@(c){(0:length(c)-1)/obj.env_Fs},envelopes);
                windows = cellfun(@(c){min(c/ec,1).*min(fliplr(c)/ec,1)},times);
                envelopes = cellfun(@(c,d){c.*d},envelopes,windows);
            end
        end

        function [spectra,freqs] = extractSpectra(obj,envelopes)
            %get signal lengths
            envLengths = cellfun('length',envelopes);

            envelopes = cellfun(@(c){c-mean(c)},envelopes);
            envelopes = cellfun(@(c){c/max(abs(c))},envelopes);
            
            envelopes = obj.attenuateEdges(envelopes);

            N = obj.spec_Nfft;

            %moving average spectral smoothing filter size
            L = fix(N*obj.spec_SmoothBw/obj.env_Fs);

            if any(envLengths>N)
                warning('extractSpectra:fftUndersampled',...
                    ['envelope length exceeds number of spectral coefficients.\n' ...
                    'Spectral metrics may be unreliable.\n'...
                    'Increase nfft, decrease signal lengths, or increase downsampling factor.\n']);                
            end
            envelopesPadded = cellfun(@(c,d){[c(:)' zeros(1,N - d)]},envelopes,num2cell(envLengths));
            
            spectra = cellfun(@(c){(abs(fft(c,N)).^2)/N},envelopesPadded);
            spectra = cellfun(@(c){2*(c(1:N/2))},spectra);
            freqs = obj.env_Fs*(0:N/2-1)/N;

            %treat spectrum as periodic for smoothing
            spectra = cellfun(@(c){[fliplr(c) c fliplr(c)]},spectra);
            spectra = cellfun(@(c){filter((1/L)*ones(1,L),1,c)},spectra);

            spectra = cellfun(@(c){c((N/2)+1:(N/2)+(N/2))},spectra);

        end

        function [metrics] = spectralMetrics(obj,spectra,freqs)

            if nargin<3
                envelopes = obj.extractEnvelopes();
                [spectra,freqs] = obj.extractSpectra(envelopes);
            end

            for i=1:size(obj.spec_PowerBins,1)
                binIxs = (freqs>=obj.spec_PowerBins(i,1) & freqs<obj.spec_PowerBins(i,2));
                binPowers(:,i) = cellfun(@(c)sum(c(binIxs)),spectra);
            end
            for i=1:size(binPowers,2)-1
                metrics.("sbpr_" + i) = binPowers(:,i)./binPowers(:,i+1);
            end

            for i=1:size(obj.spec_CentroidBins,1)
                binIxs = (freqs>=obj.spec_CentroidBins(i,1) & freqs<obj.spec_CentroidBins(i,2));
                metrics.("scntr_" + i) = cellfun(@(c)sum(c(binIxs).*(freqs(binIxs))/sum(c(binIxs))),spectra(:));
            end
            metrics = struct2table(metrics);
        end

        function [metrics,imfs,imfw,report] = emdMetrics(obj,envelopes)

            if nargin==1
                envelopes=obj.extractEnvelopes();
            end

            envelopes = obj.attenuateEdges(envelopes);

            [imfs,imfw,report] = obj.getImfs(envelopes);

            for i=1:obj.emd_MaxImf
                metrics.("sumpow_imf"+i) = cellfun(@(c)sum(abs(c(i,:))),imfs(:));
                metrics.("pow_imf"+i) = cellfun(@(c)sum(abs(c(i,:)))*(obj.env_Fs/sum(~isnan(c(i,:)))),imfs(:));
                metrics.("mu_w"+i) = cellfun(@(c)nanmean(c(i,:)),imfw(:));
                metrics.("var_w"+i) = cellfun(@(c)nanvar(c(i,:)),imfw(:));
                metrics.("sd_w"+i) = cellfun(@(c)nanstd(c(i,:)),imfw(:));
            end

            if obj.emd_MaxImf>1
                for i=1:obj.emd_MaxImf-1
                    metrics.("imf_ratio"+(i+1)+i) = metrics.("sumpow_imf"+(i+1))./metrics.("sumpow_imf"+i);
                end            
            end

            metrics = struct2table(metrics);

        end

        function [imfs,imfw,report] = getImfs(obj,envelopes)
            
            imfs = cellfun(@(c){emd(c,"SiftRelativeTolerance", ...
                obj.emd_SiftRelTol,'MaxNumIMF',obj.emd_MaxImf)},envelopes); 
            
            numImfs = cellfun(@(c)size(c,2),imfs);

            missingImfs = zeros(1,obj.emd_MaxImf);
            for i=1:obj.emd_MaxImf
                missingImfs(i) = sum(numImfs<i);
            end

            for i=1:length(imfs)
                if numImfs(i)<obj.emd_MaxImf
                    imfs{i}(:,numImfs(i)+1:obj.emd_MaxImf) = nan;
                end
            end

            for i=1:length(imfs)
                w{i} = nan(size(imfs{i},1),obj.emd_MaxImf);
                [~,~,times{i},w{i}(:,1:numImfs(i))] = hht(imfs{i}(:,1:numImfs(i)),obj.env_Fs);
            end
            w_all = vertcat(w{:});
            nan_missing = sum(isnan(w_all));            

            %replace edge values with nan
            if ~isempty(obj.emd_EdgeNull) && ~isnan(obj.emd_EdgeNull)
                n = round(obj.env_Fs*obj.emd_EdgeNull);                
                for i=1:length(w)                    
                    w{i}([1:n end-n+1:end],:) = nan;                    
                end
            end

            w_lens = cellfun(@(c)length(c),times);
            w_all = vertcat(w{:});
            nan_edge = sum(isnan(w_all)) - nan_missing;

            %replace out-of-range frequencies with 
            nan_outofrange = 0;
            if ~isempty(obj.emd_ImfFreqBounds)
                w_all(w_all<obj.emd_ImfFreqBounds(1)) = nan;
                w_all(w_all>obj.emd_ImfFreqBounds(2)) = nan;
                nan_outofrange = sum(isnan(w_all)) - nan_edge - nan_missing;
            end

            %percentile exclusions
            nan_prctileexc = 0;
            if ~isempty(obj.emd_FreqExclusionPercentile) && ~isnan(obj.emd_FreqExclusionPercentile)
                p_imf = prctile(w_all,obj.emd_FreqExclusionPercentile);
                w_all(w_all>p_imf) = nan;
                nan_prctileexc = sum(isnan(w_all)) - nan_outofrange - nan_edge - nan_missing;
            end

            report = table("imf" + (1:obj.emd_MaxImf)', ...
                missingImfs',...
                nan_outofrange'/size(w_all,1), ...
                nan_prctileexc'/size(w_all,1),...
                'VariableNames',{'imf' 'missing' 'outOfRange' 'percentileExclusions'});

            if obj.verbose
                fprintf('\n');
                disp(report);
            end

            c=0;
            imfw = cell(length(w),1);
            for i=1:length(w)
                imfw{i} = w_all(c + (1:w_lens(i)),:);
                c = c + w_lens(i);
            end

            imfw = cellfun(@(c){c'},imfw);
            imfs = cellfun(@(c){c'},imfs);

        end
    end
end