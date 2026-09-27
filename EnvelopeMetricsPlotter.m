classdef EnvelopeMetricsPlotter

    properties (Constant)
        FigureSizeInches = [0 0 6.4 4.8] % 640x480 px at 100 dpi
        ExteriorMargins = [0.09 0.09 0.03 0.08] % left, bottom, right, top - leaves room for title
    end

    methods (Static)

        function fig = plotWaveformEnvelope(t, x, t_env, env)
            % Overlay a waveform with its extracted amplitude envelope.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            plot(ax, t, x, 'Color', [0.6 0.6 0.6]);
            plot(ax, t_env, env, 'LineWidth', 1.5);
            xlabel(ax, 'time (s)');
            ylabel(ax, 'amplitude');
            legend(ax, 'waveform', 'envelope', 'Location', 'best');
            title(ax, 'Vocalic energy amplitude envelope');
            grid(ax, 'on');
        end

        function fig = plotSpectrum(freqs, spectrum)
            % Plot the envelope power spectrum.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            plot(ax, freqs, spectrum, 'LineWidth', 1.5);
            xlim(ax, [0 15]);
            xlabel(ax, 'frequency (Hz)');
            ylabel(ax, 'power');
            title(ax, 'Envelope power spectrum');
            grid(ax, 'on');
        end

        function fig = plotImfs(t_env, env, imfs)
            % Overlay the envelope with its empirical mode decomposition components.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            plot(ax, t_env, env, 'Color', [0.6 0.6 0.6]);
            plot(ax, t_env, imfs, 'LineWidth', 1.2);
            xlabel(ax, 'time (s)');
            ylabel(ax, 'amplitude');
            legendLabels = ["envelope" "imf" + (1:size(imfs, 1))];
            legend(ax, legendLabels, 'Location', 'best');
            title(ax, 'Envelope and intrinsic mode functions');
            grid(ax, 'on');
        end

        function fig = plotImfFrequencies(t_env, env, imfw)
            % Overlay the envelope with the instantaneous frequency of each IMF.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            yyaxis(ax, 'left');
            plot(ax, t_env, env, 'Color', [0.6 0.6 0.6]);
            ylabel(ax, 'envelope amplitude');
            yyaxis(ax, 'right');
            plot(ax, t_env, imfw, 'LineWidth', 1.2);
            ylabel(ax, 'instantaneous frequency (Hz)');
            xlabel(ax, 'time (s)');
            legendLabels = ["envelope" "imf" + (1:size(imfw, 1)) + " freq."];
            legend(ax, legendLabels, 'Location', 'best');
            title(ax, 'IMF instantaneous frequency');
            grid(ax, 'on');
        end

        function saveFigure(fig, filepath)
            % Export an stFig to a PNG with a consistent, controlled resolution.
            [d, f] = fileparts(filepath);
            set(0, 'CurrentFigure', fig.Handle);
            fig.print(fullfile(d, f), 'resolution', '-r150');
        end

    end

    methods (Static, Access = private)

        function fig = newFigure()
            fig = stFig([1 1], EnvelopeMetricsPlotter.ExteriorMargins, [], 'Theme', 'light');
            set(fig.Handle, 'Units', 'inches', 'Position', EnvelopeMetricsPlotter.FigureSizeInches);
        end

    end

end
