%% EnvelopeMetrics demo
% Runnable walkthrough of the EnvelopeMetrics pipeline. Mirrors README.md and
% TECHNICAL.md, and regenerates every figure those documents embed.
%
% Run with: matlab -batch "run('demo.m')"

repoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(repoDir, 'utils'));

figDir = fullfile(repoDir, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Calculating the metrics
% Load audio files (sampling rates must match across files) and get the full
% metrics table in one call.
files = ["data/example1.wav" "data/example2.wav"];
[X, Fs] = arrayfun(@(c)audioread(c), files, 'UniformOutput', false);

em = envelopeMetrics(X, Fs{1});
metrics = em.getMetrics();
disp(metrics);

%% The vocalic energy amplitude envelope
% The envelope is a "vocalic energy amplitude envelope": the low-pass
% zero-phase-filtered absolute value of the band-pass filtered waveform. It
% fits the magnitude of the vocalic energy waveform, not the raw waveform.
x = X{1} / max(abs(X{1}));
fs = Fs{1};
t = (0:length(x) - 1) / fs;

[env, t_env] = em.extractEnvelopes();

fig1 = EnvelopeMetricsPlotter.plotWaveformEnvelope(t, x, t_env{1}, env{1});
EnvelopeMetricsPlotter.saveFigure(fig1, fullfile(figDir, 'waveform_envelope.png'));

%% Envelope spectrum metrics
% The envelope is zero-centered, edge-attenuated, and zero-padded before FFT.
[spectra, freqs] = em.extractSpectra(env);

fig2 = EnvelopeMetricsPlotter.plotSpectrum(freqs, spectra{1});
EnvelopeMetricsPlotter.saveFigure(fig2, fullfile(figDir, 'spectrum.png'));

psMetrics = em.spectralMetrics(spectra, freqs);
disp(psMetrics);

%% Empirical mode decomposition metrics
% The first IMF captures syllable-timescale oscillations, the second
% stress-timescale oscillations. Higher-order IMFs capture lower-frequency,
% phrase-timescale oscillations for longer chunks.
[imfs, imfw] = em.getImfs(env);

fig3 = EnvelopeMetricsPlotter.plotImfs(t_env{1}, env{1}, imfs{1}(1:2, :));
EnvelopeMetricsPlotter.saveFigure(fig3, fullfile(figDir, 'imfs.png'));

fig4 = EnvelopeMetricsPlotter.plotImfFrequencies(t_env{1}, imfw{1}(1:2, :));
EnvelopeMetricsPlotter.saveFigure(fig4, fullfile(figDir, 'imf_freq.png'));

% Combined overview: the envelope as fed into EMD (attenuated), its IMFs,
% and their instantaneous frequencies, stacked and x-aligned in one figure.
attenuatedEnv = em.attenuateEdges(env);
fig5 = EnvelopeMetricsPlotter.plotEmdOverview(t_env{1}, attenuatedEnv{1}, imfs{1}(1:2, :), imfw{1}(1:2, :));
EnvelopeMetricsPlotter.saveFigure(fig5, fullfile(figDir, 'emd_overview.png'));

emdMetrics = em.emdMetrics(env);
disp(emdMetrics);

close([fig1.Handle fig2.Handle fig3.Handle fig4.Handle fig5.Handle]);

%% Parameters and presets
% getParams/setParams store/update every tunable property in a plain struct,
% so a corpus-specific configuration can be saved and reapplied later.
preset = envelopeMetrics.defaultParams();
preset.env_Passband = [300 3400];
preset.emd_MaxImf = 4;

em2 = envelopeMetrics(X, Fs{1}).setParams(preset);
fprintf('em2 uses %d IMFs and a %g-%g Hz passband\n', ...
    em2.emd_MaxImf, em2.env_Passband(1), em2.env_Passband(2));

%% Table input
% A table of file paths (and optional per-token t0/t1 start/end times, in
% seconds) can be used instead of raw audio. Missing files are skipped with
% a warning. getMetrics() then returns the input table with an added Fs
% column, each token's envelope, and every metric as trailing columns.
T = readtable(fullfile(repoDir, 'data', 'exampleTokensTable.csv'));
em3 = envelopeMetrics(T);
resultsTable = em3.getMetrics();
disp(resultsTable);
