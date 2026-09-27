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
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            plot(ax, t, x, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            plot(ax, t_env, env, 'LineWidth', 2.5);
            xlabel(ax, 'time (s)', 'FontSize', 14);
            ylabel(ax, 'amplitude', 'FontSize', 14);
            legend(ax, 'waveform', 'envelope', 'Location', 'best', 'FontSize', 12);
            title(ax, 'Vocalic energy amplitude envelope', 'FontSize', 15);
            grid(ax, 'on');
        end

        function fig = plotSpectrum(freqs, spectrum)
            % Plot the envelope power spectrum.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            plot(ax, freqs, spectrum, 'LineWidth', 2.5);
            xlim(ax, [0 15]);
            xlabel(ax, 'frequency (Hz)', 'FontSize', 14);
            ylabel(ax, 'power', 'FontSize', 14);
            title(ax, 'Envelope power spectrum', 'FontSize', 15);
            grid(ax, 'on');
        end

        function fig = plotImfs(t_env, env, imfs)
            % Overlay the envelope with its empirical mode decomposition components.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            plot(ax, t_env, env, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            plot(ax, t_env, imfs, 'LineWidth', 2);
            xlabel(ax, 'time (s)', 'FontSize', 14);
            ylabel(ax, 'amplitude', 'FontSize', 14);
            legendLabels = ["envelope" "imf" + (1:size(imfs, 1))];
            legend(ax, legendLabels, 'Location', 'best', 'FontSize', 12);
            title(ax, 'Envelope and intrinsic mode functions', 'FontSize', 15);
            grid(ax, 'on');
        end

        function fig = plotImfFrequencies(t_env, env, imfw)
            % Overlay the envelope with the instantaneous frequency of each IMF.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            yyaxis(ax, 'left');
            plot(ax, t_env, env, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            ylabel(ax, 'envelope amplitude', 'FontSize', 14);
            yyaxis(ax, 'right');
            plot(ax, t_env, imfw, 'LineWidth', 2);
            ylabel(ax, 'instantaneous frequency (Hz)', 'FontSize', 14);
            xlabel(ax, 'time (s)', 'FontSize', 14);
            legendLabels = ["envelope" "imf" + (1:size(imfw, 1)) + " freq."];
            legend(ax, legendLabels, 'Location', 'best', 'FontSize', 12);
            title(ax, 'IMF instantaneous frequency', 'FontSize', 15);
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

        function setAxisFontSize(ax)
            % Tick label size; title/axis/legend labels set their own FontSize above.
            set(ax, 'FontSize', 13);
        end

    end

end
