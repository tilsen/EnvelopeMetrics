classdef EnvelopeMetricsPlotter

    properties (Constant)
        FigureSizeInches = [0 0 6.4 4.8] % 640x480 px at 100 dpi
        ExteriorMargins = [0.22 0.22 0.09 0.17] % left, bottom, right, top - room for larger fonts
    end

    methods (Static)

        function fig = plotWaveformEnvelope(t, x, t_env, env)
            % Overlay a waveform with its extracted amplitude envelope.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            plot(ax, t, x, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            plot(ax, t_env, env, 'LineWidth', 2);
            xlabel(ax, 'time (s)', 'FontSize', 20);
            ylabel(ax, 'amplitude', 'FontSize', 20);
            legend(ax, 'waveform', 'envelope', 'Location', 'best', 'FontSize', 15);
            title(ax, 'Vocalic energy amplitude envelope', 'FontSize', 22);
            grid(ax, 'on');
        end

        function fig = plotSpectrum(freqs, spectrum)
            % Plot the envelope power spectrum.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            plot(ax, freqs, spectrum, 'LineWidth', 2);
            xlim(ax, [0 15]);
            xlabel(ax, 'frequency (Hz)', 'FontSize', 20);
            ylabel(ax, 'power', 'FontSize', 20);
            title(ax, 'Envelope power spectrum', 'FontSize', 22);
            grid(ax, 'on');
        end

        function fig = plotImfs(t_env, env, imfs)
            % Overlay the envelope with its empirical mode decomposition components.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            plot(ax, t_env, env, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            plot(ax, t_env, imfs, 'LineWidth', 2);
            xlabel(ax, 'time (s)', 'FontSize', 20);
            ylabel(ax, 'amplitude', 'FontSize', 20);
            legendLabels = ["envelope" "imf" + (1:size(imfs, 1))];
            legend(ax, legendLabels, 'Location', 'best', 'FontSize', 15);
            title(ax, 'Envelope and intrinsic mode functions', 'FontSize', 22);
            grid(ax, 'on');
        end

        function fig = plotImfFrequencies(t_env, imfw)
            % Plot the instantaneous frequency of each IMF.
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax);
            plot(ax, t_env, imfw, 'LineWidth', 2);
            xlabel(ax, 'time (s)', 'FontSize', 20);
            ylabel(ax, 'instantaneous frequency (Hz)', 'FontSize', 20);
            legendLabels = "imf" + (1:size(imfw, 1)) + " freq.";
            legend(ax, legendLabels, 'Location', 'best', 'FontSize', 15);
            title(ax, 'IMF instantaneous frequency', 'FontSize', 22);
            grid(ax, 'on');
        end

        function fig = plotEmdOverview(t_env, envelope, imfs, imfw)
            % Vertically stacked, x-aligned panels: the envelope as fed into
            % EMD (after padding/recentering/attenuation), its intrinsic
            % mode functions, and their instantaneous frequencies.
            fig = stFig([3 1], EnvelopeMetricsPlotter.ExteriorMargins, [0 0.01], Theme='light', Aspect=0.85);
            set(fig.Handle, 'Units', 'inches');

            xl = [min(t_env) max(t_env)];

            ax1 = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax1);
            plot(ax1, t_env, envelope, 'LineWidth', 2);
            ylabel(ax1, 'amplitude', 'FontSize', 20);
            title(ax1, 'Envelope, intrinsic mode functions, and instantaneous frequency', 'FontSize', 22);
            xlim(ax1, xl);
            set(ax1, 'XTickLabel', []);
            grid(ax1, 'on');

            ax2 = fig.Axes(2);
            EnvelopeMetricsPlotter.setAxisFontSize(ax2);
            plot(ax2, t_env, imfs, 'LineWidth', 2);
            ylabel(ax2, 'amplitude', 'FontSize', 20);
            legendLabels = "imf" + (1:size(imfs, 1));
            legend(ax2, legendLabels, 'Location', 'best', 'FontSize', 15);
            xlim(ax2, xl);
            set(ax2, 'XTickLabel', []);
            grid(ax2, 'on');

            ax3 = fig.Axes(3);
            EnvelopeMetricsPlotter.setAxisFontSize(ax3);
            plot(ax3, t_env, imfw, 'LineWidth', 2);
            xlabel(ax3, 'time (s)', 'FontSize', 20);
            ylabel(ax3, 'inst. freq. (Hz)', 'FontSize', 20);
            legendLabels = "imf" + (1:size(imfw, 1)) + " freq.";
            legend(ax3, legendLabels, 'Location', 'best', 'FontSize', 15);
            xlim(ax3, xl);
            grid(ax3, 'on');
        end

        function saveFigure(fig, filepath)
            % Export to a PNG via exportgraphics rather than stFig's own
            % print()/printpng (a thin wrapper over the legacy print()
            % hardcopy path).
            exportgraphics(fig.Handle, filepath, 'Resolution', 150);
        end

    end

    methods (Static, Access = private)

        function fig = newFigure()
            fig = stFig([1 1], EnvelopeMetricsPlotter.ExteriorMargins, [], 'Theme', 'light');
            set(fig.Handle, 'Units', 'inches', 'Position', EnvelopeMetricsPlotter.FigureSizeInches);
        end

        function setAxisFontSize(ax)
            % Tick label size; title/axis/legend labels set their own FontSize above.
            set(ax, 'FontSize', 18);
        end

    end

end
