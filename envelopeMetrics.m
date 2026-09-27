classdef envelopeMetrics

    properties
        data % input table (file/t0/t1), when constructed from one; empty otherwise
        X % per-token waveforms, one row vector per cell (valid tokens only)
        Fs % common sampling rate of X, in Hz
    end

    properties (Access = private)
        valid_ = [] % logical, one per row of data; false = file/channel skipped (getMetrics() fills that row with nan)
    end

    properties (Dependent)
        env_Fs % = Fs / env_Downsample; recomputed on access, never stale
    end

    properties
        env_Passband (1,2) double {mustBeValidPassband} = [400 4000] % [low high]; high may be nan to mean Nyquist
        env_Lowpass (1,1) double {mustBePositive} = 10
        env_BandpassFilterOrder (1,1) double {mustBePositive, mustBeInteger} = 4
        env_LowpassFilterOrder (1,1) double {mustBePositive, mustBeInteger} = 4
        env_Downsample (1,1) double {mustBePositive, mustBeInteger} = 100
        env_Rescale (1,1) logical = true
        env_ZeroPad (1,1) double {mustBeNonnegative} = 0 % seconds of zero-padding appended to each end of the waveform before filtering, then trimmed back off after the low-pass filter; reduces filtfilt edge transients for short tokens
        env_TukeywinParam (1,1) double {mustBeNonnegativeOrNan} = nan
        env_EdgeAttenutation (1,1) double {mustBeNonnegativeOrNan} = 0.05 %period of time to attenuate edges
        spec_Nfft (1,1) double {mustBePositive, mustBeInteger} = 2048
        spec_SmoothBw (1,1) double {mustBePositive} = 1
        spec_PowerBins (:,2) double {mustBeValidBinRows} = [1 3.5; 3.5 10]
        spec_CentroidBins (:,2) double {mustBeValidBinRows} = [1 10]
        emd_MaxImf (1,1) double {mustBePositive, mustBeInteger} = 3
        emd_EdgeNull (1,1) double {mustBeNonnegativeOrNan} = 0.1
        emd_SiftRelTol (1,1) double {mustBePositive} = 0.1
        emd_SiftMaxIterations (1,1) double {mustBeNonnegativeOrNan} = nan % nan = use emd()'s default (100)
        emd_MaxNumExtrema (1,1) double {mustBeNonnegativeOrNan} = nan % nan = use emd()'s default (1)
        emd_MaxEnergyRatio (1,1) double {mustBeNonnegativeOrNan} = nan % nan = use emd()'s default (20)
        emd_Interpolation (1,1) string = "" % "" = use emd()'s default ('spline'); or "pchip"
        emd_HhtFrequencyLimits (1,:) double {mustBeEmptyOrIncreasingPair} = [] % [] = use hht()'s default ([0 env_Fs/2])
        emd_FreqExclusionPercentile (1,1) double {mustBeInRange0to100OrNan} = 99
        emd_AutoImfFreqBoundsDb (1,1) double {mustBeNegative} = -10 % dB point of the envelope lowpass filter used to auto-derive emd_ImfFreqBounds
        includeEnvelope (1,1) logical = false % if true, getMetrics() appends envFs and envelope as the last two columns of its output table
        verbose (1,1) logical = false
    end

    properties (Access = private)
        emd_ImfFreqBounds_ = [] % backing store; [] = auto-derive from env_Lowpass/env_LowpassFilterOrder/emd_AutoImfFreqBoundsDb
    end

    properties (Dependent)
        emd_ImfFreqBounds % [low high]; [] (default) auto-derives from the envelope lowpass filter - see emd_AutoImfFreqBoundsDb
    end

    methods
        function obj = envelopeMetrics(data,Fs)
            arguments
                data = []
                Fs = []
            end
            if nargin==0, return; end

            if istable(data)
                obj = obj.initFromTable(data);
            else
                obj = obj.initFromAudio(data,Fs);
            end
        end

        function fs = get.env_Fs(obj)
            fs = obj.Fs/obj.env_Downsample;
        end

        function bounds = get.emd_ImfFreqBounds(obj)
            if isempty(obj.emd_ImfFreqBounds_)
                bounds = obj.autoImfFreqBounds();
            else
                bounds = obj.emd_ImfFreqBounds_;
            end
        end

        function obj = set.emd_ImfFreqBounds(obj,bounds)
            if ~isempty(bounds)
                mustBeEmptyOrIncreasingPair(bounds);
            end
            obj.emd_ImfFreqBounds_ = bounds;
        end

        function bounds = autoImfFreqBounds(obj)
            % The frequency at which the envelope's lowpass filter (env_Lowpass,
            % env_LowpassFilterOrder) reaches emd_AutoImfFreqBoundsDb, using the
            % closed-form Butterworth magnitude response (independent of Fs).
            n = obj.env_LowpassFilterOrder;
            dB = obj.emd_AutoImfFreqBoundsDb;
            ratio = (10^(-dB/10) - 1)^(1/(2*n));
            bounds = [0 obj.env_Lowpass*ratio];
        end

        function [T,envelopes,times,spectra,freqs,imfs,imfw] = getMetrics(obj)
            % A single output table: starts from the input table (file, t0,
            % t1) when this object was constructed from one, otherwise one
            % row per token; adds a common Fs column and appends every
            % metric as trailing columns. If includeEnvelope is true, envFs
            % and each token's envelope are appended last. Rows whose file
            % or channel was invalid (see initFromTable) are retained, with
            % nan for Fs, every metric, and (if included) envFs/envelope.
            [envelopes,times,~,~,obj] = obj.extractEnvelopes();
            [spectra, freqs] = obj.extractSpectra(envelopes);
            spectralMetrics = obj.spectralMetrics(spectra,freqs);
            [emdMetrics,imfs,imfw] = obj.emdMetrics(envelopes);
            metrics = [spectralMetrics emdMetrics];
            cols = metrics.Properties.VariableNames;
            colsOrdered = cols(~cellfun('isempty',(regexp(cols,'sbpr|scntr|imf_ratio', 'once'))));
            colsOrdered = [colsOrdered sort(setdiff(cols,colsOrdered))];
            metrics = metrics(:,colsOrdered);

            if istable(obj.data) && width(obj.data)>0
                T = obj.data;
                valid = obj.valid_;
            else
                T = table((1:numel(obj.X))','VariableNames',{'token'});
                valid = true(height(T),1);
            end
            n = height(T);

            Fs_col = nan(n,1);
            Fs_col(valid) = obj.Fs;
            T.Fs = Fs_col;

            metricsFull = array2table(nan(n,width(metrics)),'VariableNames',metrics.Properties.VariableNames);
            metricsFull(valid,:) = metrics;
            T = [T metricsFull];

            if obj.includeEnvelope
                envFs_col = nan(n,1);
                envFs_col(valid) = obj.env_Fs;
                T.envFs = envFs_col;

                envelopeCol = cell(n,1);
                envelopeCol(valid) = envelopes(:);
                T.envelope = envelopeCol;
            end
        end

        function [envelopes, times,passbandFiltered, lowpassFiltered,obj] = extractEnvelopes(obj)
            x = cellfun(@(c){c-mean(c)},obj.X);

            padN = round(obj.env_ZeroPad*obj.Fs);
            if padN>0
                x = cellfun(@(c){[zeros(1,padN) c zeros(1,padN)]},x);
            end

            passband = obj.env_Passband;
            if isnan(passband(2))
                passband(2) = obj.Fs/2;
            end

            [bbp,abp] = butter(obj.env_BandpassFilterOrder,passband/(obj.Fs/2));
            [blp,alp] = butter(obj.env_LowpassFilterOrder,obj.env_Lowpass/(obj.Fs/2));

            passbandFiltered = cellfun(@(c){filtfilt(bbp,abp,c)},x);
            lowpassFiltered = cellfun(@(c){filtfilt(blp,alp,abs(c))},passbandFiltered);

            if padN>0
                lowpassFiltered = cellfun(@(c){c(padN+1:end-padN)},lowpassFiltered);
            end

            envelopes = cellfun(@(c){downsample(c,obj.env_Downsample)},lowpassFiltered);

            if obj.env_Rescale
                envelopes = cellfun(@(c){c/max(abs(c))},envelopes);
            end
            times = cellfun(@(c){(0:length(c)-1)/obj.env_Fs},envelopes);

        end

        function [envelopes] = attenuateEdges(obj,envelopes)
            if ~isnan(obj.env_TukeywinParam)
                envelopes = cellfun(@(c){c.*tukeywin(length(c),obj.env_TukeywinParam)'},envelopes);
            end
            if ~isnan(obj.env_EdgeAttenutation)
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

            if size(obj.spec_PowerBins,1)>1 && ...
                    any(obj.spec_PowerBins(2:end,1) ~= obj.spec_PowerBins(1:end-1,2))
                warning('envelopeMetrics:NonContiguousPowerBins', ...
                    ['spec_PowerBins rows are not contiguous (row i''s high edge should ' ...
                    'equal row i+1''s low edge); sbpr_i (power ratio of bin i to bin i+1) ' ...
                    'may not be a meaningful adjacent-band comparison.']);
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

            emdArgs = obj.emdOptionArgs();
            imfs = cellfun(@(c){emd(c,emdArgs{:})},envelopes);

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

            hhtArgs = obj.hhtOptionArgs();
            for i=1:length(imfs)
                w{i} = nan(size(imfs{i},1),obj.emd_MaxImf);
                [~,~,times{i},w{i}(:,1:numImfs(i))] = hht(imfs{i}(:,1:numImfs(i)),obj.env_Fs,hhtArgs{:});
            end
            w_all = vertcat(w{:});
            nan_missing = sum(isnan(w_all));

            %replace edge values with nan
            if ~isnan(obj.emd_EdgeNull)
                n = round(obj.env_Fs*obj.emd_EdgeNull);
                for i=1:length(w)
                    w{i}([1:n end-n+1:end],:) = nan;
                end
            end

            w_lens = cellfun(@(c)length(c),times);
            w_all = vertcat(w{:});
            nan_edge = sum(isnan(w_all)) - nan_missing;

            %replace out-of-range frequencies with nan
            bounds = obj.emd_ImfFreqBounds;
            w_all(w_all<bounds(1)) = nan;
            w_all(w_all>bounds(2)) = nan;
            nan_outofrange = sum(isnan(w_all)) - nan_edge - nan_missing;

            %percentile exclusions
            nan_prctileexc = 0;
            if ~isnan(obj.emd_FreqExclusionPercentile)
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

        function params = getParams(obj)
            % Every tunable analysis parameter as a plain struct, for saving
            % alongside results (provenance) or as a reusable corpus-specific
            % preset (see setParams). emd_ImfFreqBounds is captured as its raw
            % setting (possibly [] for "auto"), not the resolved value, so a
            % preset built this way keeps auto-deriving it when applied to an
            % object with different envelope-filter properties.
            names = envelopeMetrics.paramNames();
            params = struct();
            for i = 1:numel(names)
                name = names{i};
                if strcmp(name,'emd_ImfFreqBounds')
                    params.(name) = obj.emd_ImfFreqBounds_;
                else
                    params.(name) = obj.(name);
                end
            end
        end

        function obj = setParams(obj,params)
            % Apply a struct of tunable parameters (as returned by getParams,
            % or loaded from a saved corpus-specific preset) to this object.
            fields = fieldnames(params);
            valid = envelopeMetrics.paramNames();
            for i = 1:numel(fields)
                if ~ismember(fields{i},valid)
                    error('envelopeMetrics:InvalidParam','Unknown parameter: %s',fields{i});
                end
                obj.(fields{i}) = params.(fields{i});
            end
        end

    end

    methods (Access = private)

        function obj = initFromAudio(obj,X,Fs)
            if ~iscell(X) || ~all(cellfun(@(c) isnumeric(c) && isvector(c),X))
                error('envelopeMetrics:InvalidInput', ...
                    'Input waveforms must be provided as a cell array of numeric row or column vectors.');
            end
            if isempty(Fs) || ~isscalar(Fs) || ~isnumeric(Fs) || Fs<=0
                error('envelopeMetrics:InvalidInput', ...
                    'Sampling rate must be provided as a positive scalar second input.');
            end
            obj.X = cellfun(@(c){c(:)'},X);
            obj.Fs = Fs;
            obj.data = table();
            obj.valid_ = true(numel(X),1);
        end

        function obj = initFromTable(obj,T)
            % Table input: one row per token, with a 'file' column of audio
            % paths and optional 't0'/'t1' columns (start/end time in
            % seconds; default to the whole file when absent or nan) and an
            % optional 'channel' column (1-based; defaults to 1). Rows whose
            % file doesn't exist or whose channel is invalid are kept in
            % data (with a warning), but excluded from X; getMetrics() fills
            % those rows with nan rather than dropping them.
            if ~ismember('file',T.Properties.VariableNames)
                error('envelopeMetrics:InvalidInput', ...
                    'Table input must contain a ''file'' column of audio file paths.');
            end
            if ~ismember('t0',T.Properties.VariableNames)
                T.t0 = nan(height(T),1);
            end
            if ~ismember('t1',T.Properties.VariableNames)
                T.t1 = nan(height(T),1);
            end
            if ~ismember('channel',T.Properties.VariableNames)
                T.channel = ones(height(T),1);
            else
                T.channel(isnan(T.channel)) = 1;
            end

            files = string(T.file);
            exists = isfile(files);
            for i = find(~exists)'
                warning('envelopeMetrics:FileNotFound', ...
                    'File does not exist and will be skipped: %s',files(i));
            end

            valid = exists;
            for i = find(exists)'
                numChannels = audioinfo(files(i)).NumChannels;
                if T.channel(i)<1 || T.channel(i)>numChannels
                    valid(i) = false;
                    warning('envelopeMetrics:InvalidChannel', ...
                        'File has %d channel(s); requested channel %d does not exist and will be skipped: %s', ...
                        numChannels,T.channel(i),files(i));
                end
            end

            if ~any(valid)
                error('envelopeMetrics:InvalidInput', ...
                    'None of the files in the input table exist with a valid channel.');
            end

            validIx = find(valid);
            X = cell(numel(validIx),1);
            fs = nan(numel(validIx),1);
            for k = 1:numel(validIx)
                i = validIx(k);
                [y,fsi] = audioread(files(i));
                y = y(:,T.channel(i));
                fs(k) = fsi;

                i0 = 1;
                i1 = size(y,1);
                if ~isnan(T.t0(i))
                    i0 = max(1,round(T.t0(i)*fsi)+1);
                end
                if ~isnan(T.t1(i))
                    i1 = min(size(y,1),round(T.t1(i)*fsi));
                end
                X{k} = y(i0:i1)';
            end

            if any(fs ~= fs(1))
                error('envelopeMetrics:InvalidInput', ...
                    'Sampling rates must be identical across all files in the input table.');
            end

            obj.X = X;
            obj.Fs = fs(1);
            obj.data = T;
            obj.valid_ = valid;
        end

        function args = emdOptionArgs(obj)
            args = {'SiftRelativeTolerance', obj.emd_SiftRelTol, 'MaxNumIMF', obj.emd_MaxImf};
            if ~isnan(obj.emd_SiftMaxIterations)
                args = [args, {'SiftMaxIterations', obj.emd_SiftMaxIterations}];
            end
            if ~isnan(obj.emd_MaxNumExtrema)
                args = [args, {'MaxNumExtrema', obj.emd_MaxNumExtrema}];
            end
            if ~isnan(obj.emd_MaxEnergyRatio)
                args = [args, {'MaxEnergyRatio', obj.emd_MaxEnergyRatio}];
            end
            if strlength(obj.emd_Interpolation) > 0
                args = [args, {'Interpolation', obj.emd_Interpolation}];
            end
        end

        function args = hhtOptionArgs(obj)
            args = {};
            if ~isempty(obj.emd_HhtFrequencyLimits)
                args = {'FrequencyLimits', obj.emd_HhtFrequencyLimits};
            end
        end

    end

    methods (Static)

        function names = paramNames()
            names = {'env_Passband','env_Lowpass','env_BandpassFilterOrder', ...
                'env_LowpassFilterOrder','env_Downsample','env_Rescale','env_ZeroPad', ...
                'env_TukeywinParam','env_EdgeAttenutation', ...
                'spec_Nfft','spec_SmoothBw','spec_PowerBins','spec_CentroidBins', ...
                'emd_MaxImf','emd_EdgeNull','emd_SiftRelTol', ...
                'emd_SiftMaxIterations','emd_MaxNumExtrema','emd_MaxEnergyRatio', ...
                'emd_Interpolation','emd_HhtFrequencyLimits', ...
                'emd_ImfFreqBounds','emd_AutoImfFreqBoundsDb', ...
                'emd_FreqExclusionPercentile','includeEnvelope','verbose'};
        end

        function params = defaultParams()
            % Convenience preset equal to a freshly constructed object's parameters.
            em = envelopeMetrics();
            params = em.getParams();
        end

    end
end

function mustBeValidPassband(v)
if numel(v)~=2
    error('envelopeMetrics:InvalidPassband','env_Passband must have exactly 2 elements: [low high].');
end
if ~isfinite(v(1)) || v(1)<0
    error('envelopeMetrics:InvalidPassband', ...
        'env_Passband(1) (the low cutoff) must be a finite, non-negative number.');
end
if ~isnan(v(2)) && v(2)<=v(1)
    error('envelopeMetrics:InvalidPassband', ...
        'env_Passband(2) (the high cutoff) must be greater than env_Passband(1), or nan to use the Nyquist frequency.');
end
end

function mustBeNonnegativeOrNan(v)
if ~isnan(v) && v<0
    error('envelopeMetrics:InvalidParam','Value must be non-negative, or nan to disable.');
end
end

function mustBeValidBinRows(v)
if isempty(v), return; end
if any(v(:,2)<=v(:,1))
    error('envelopeMetrics:InvalidBinRows','Each bin row must be [low high] with high > low.');
end
end

function mustBeInRange0to100OrNan(v)
if ~isnan(v) && (v<0 || v>100)
    error('envelopeMetrics:InvalidParam','Value must be between 0 and 100, or nan to disable.');
end
end

function mustBeEmptyOrIncreasingPair(v)
if isempty(v), return; end
if numel(v)~=2 || v(2)<=v(1)
    error('envelopeMetrics:InvalidParam','Value must be empty, or a 2-element vector [low high] with high > low.');
end
end
