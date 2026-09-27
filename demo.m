%% EnvelopeMetrics demo
% Runnable walkthrough of the EnvelopeMetrics pipeline. Mirrors README.md and
% TECHNICAL.md, and regenerates every figure those documents embed.
%
% Run with: matlab -batch "run('demo.m')"

repoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(repoDir, 'figtools'));

figDir = fullfile(repoDir, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Calculating the metrics
% Load audio files (sampling rates must match across files) and get the full
% metrics table in one call.
files = ["example.wav" "example2.wav"];
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

fig4 = EnvelopeMetricsPlotter.plotImfFrequencies(t_env{1}, env{1}, imfw{1}(1:2, :));
EnvelopeMetricsPlotter.saveFigure(fig4, fullfile(figDir, 'imf_freq.png'));

emdMetrics = em.emdMetrics(env);
disp(emdMetrics);

close([fig1.Handle fig2.Handle fig3.Handle fig4.Handle]);

%% Parameters and presets
% Every tunable property round-trips through getParams/setParams as a plain
% struct, so a corpus-specific configuration can be saved and reapplied later.
preset = envelopeMetrics.defaultParams();
preset.env_Passband = [300 3400];
preset.emd_MaxImf = 4;

em2 = envelopeMetrics(X, Fs{1}).setParams(preset);
fprintf('em2 uses %d IMFs and a %g-%g Hz passband\n', ...
    em2.emd_MaxImf, em2.env_Passband(1), em2.env_Passband(2));
