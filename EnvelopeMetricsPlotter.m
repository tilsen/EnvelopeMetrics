classdef EnvelopeMetricsPlotter

    properties (Constant)
        FigureSizeInches = [0 0 6.4 5.4] % taller than 4:3 to leave headroom for title+legend at these font sizes
        EmdOverviewSizeInches = [0 0 6.4 8.6] % same width as FigureSizeInches, taller for 3 stacked panels
        ExteriorMargins = [0.22 0.22 0.09 0.24] % left, bottom, right, top - room for larger fonts
        AxesFontSize = 12 % tick label size
        LabelFontSize = 14 % x/y axis label size
        TitleFontSize = 16
        LegendFontSize = 13
    end

    methods (Static)

        function fig = plotWaveformEnvelope(t, x, t_env, env, opts)
            % Overlay a waveform with its extracted amplitude envelope.
            arguments
                t, x, t_env, env
                opts.AxesFontSize (1,1) double = EnvelopeMetricsPlotter.AxesFontSize
                opts.LabelFontSize (1,1) double = EnvelopeMetricsPlotter.LabelFontSize
                opts.TitleFontSize (1,1) double = EnvelopeMetricsPlotter.TitleFontSize
                opts.LegendFontSize (1,1) double = EnvelopeMetricsPlotter.LegendFontSize
            end
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax, opts.AxesFontSize);
            plot(ax, t, x, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            plot(ax, t_env, env, 'LineWidth', 2);
            xlabel(ax, 'time (s)', 'FontSize', opts.LabelFontSize);
            ylabel(ax, 'amplitude', 'FontSize', opts.LabelFontSize);
            legend(ax, 'waveform', 'envelope', 'Location', 'northeast', 'FontSize', opts.LegendFontSize);
            title(ax, 'Vocalic energy amplitude envelope', 'FontSize', opts.TitleFontSize);
            grid(ax, 'on');
        end

        function fig = plotSpectrum(freqs, spectrum, opts)
            % Plot the envelope power spectrum.
            arguments
                freqs, spectrum
                opts.AxesFontSize (1,1) double = EnvelopeMetricsPlotter.AxesFontSize
                opts.LabelFontSize (1,1) double = EnvelopeMetricsPlotter.LabelFontSize
                opts.TitleFontSize (1,1) double = EnvelopeMetricsPlotter.TitleFontSize
                opts.LegendFontSize (1,1) double = EnvelopeMetricsPlotter.LegendFontSize
            end
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax, opts.AxesFontSize);
            plot(ax, freqs, spectrum, 'LineWidth', 2);
            xlim(ax, [0 15]);
            xlabel(ax, 'frequency (Hz)', 'FontSize', opts.LabelFontSize);
            ylabel(ax, 'power', 'FontSize', opts.LabelFontSize);
            title(ax, 'Envelope power spectrum', 'FontSize', opts.TitleFontSize);
            grid(ax, 'on');
        end

        function fig = plotImfs(t_env, env, imfs, opts)
            % Overlay the envelope with its empirical mode decomposition components.
            arguments
                t_env, env, imfs
                opts.AxesFontSize (1,1) double = EnvelopeMetricsPlotter.AxesFontSize
                opts.LabelFontSize (1,1) double = EnvelopeMetricsPlotter.LabelFontSize
                opts.TitleFontSize (1,1) double = EnvelopeMetricsPlotter.TitleFontSize
                opts.LegendFontSize (1,1) double = EnvelopeMetricsPlotter.LegendFontSize
            end
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax, opts.AxesFontSize);
            plot(ax, t_env, env, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            plot(ax, t_env, imfs, 'LineWidth', 2);
            xlabel(ax, 'time (s)', 'FontSize', opts.LabelFontSize);
            ylabel(ax, 'amplitude', 'FontSize', opts.LabelFontSize);
            legendLabels = ["envelope" "imf" + (1:size(imfs, 1))];
            legend(ax, legendLabels, 'Location', 'northeast', 'FontSize', opts.LegendFontSize);
            title(ax, 'Envelope and intrinsic mode functions', 'FontSize', opts.TitleFontSize);
            grid(ax, 'on');
        end

        function fig = plotImfFrequencies(t_env, imfw, opts)
            % Plot the instantaneous frequency of each IMF.
            arguments
                t_env, imfw
                opts.AxesFontSize (1,1) double = EnvelopeMetricsPlotter.AxesFontSize
                opts.LabelFontSize (1,1) double = EnvelopeMetricsPlotter.LabelFontSize
                opts.TitleFontSize (1,1) double = EnvelopeMetricsPlotter.TitleFontSize
                opts.LegendFontSize (1,1) double = EnvelopeMetricsPlotter.LegendFontSize
            end
            fig = EnvelopeMetricsPlotter.newFigure();
            ax = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax, opts.AxesFontSize);
            plot(ax, t_env, imfw, 'LineWidth', 2);
            xlabel(ax, 'time (s)', 'FontSize', opts.LabelFontSize);
            ylabel(ax, 'instantaneous frequency (Hz)', 'FontSize', opts.LabelFontSize);
            legendLabels = "imf" + (1:size(imfw, 1)) + " freq.";
            legend(ax, legendLabels, 'Location', 'northeast', 'FontSize', opts.LegendFontSize);
            title(ax, 'IMF instantaneous frequency', 'FontSize', opts.TitleFontSize);
            grid(ax, 'on');
        end

        function fig = plotEmdOverview(t_env, envelope, imfs, imfw, opts)
            % Vertically stacked, x-aligned panels: the envelope as fed into
            % EMD (after padding/recentering/attenuation), its intrinsic
            % mode functions, and their instantaneous frequencies.
            arguments
                t_env, envelope, imfs, imfw
                opts.AxesFontSize (1,1) double = EnvelopeMetricsPlotter.AxesFontSize
                opts.LabelFontSize (1,1) double = EnvelopeMetricsPlotter.LabelFontSize
                opts.TitleFontSize (1,1) double = EnvelopeMetricsPlotter.TitleFontSize
                opts.LegendFontSize (1,1) double = EnvelopeMetricsPlotter.LegendFontSize
            end
            fig = EnvelopeMetricsPlotter.createSizedFigure([3 1], EnvelopeMetricsPlotter.ExteriorMargins, ...
                [0 0.01], EnvelopeMetricsPlotter.EmdOverviewSizeInches);

            xl = [min(t_env) max(t_env)];

            ax1 = fig.Axes(1);
            EnvelopeMetricsPlotter.setAxisFontSize(ax1, opts.AxesFontSize);
            plot(ax1, t_env, envelope, 'LineWidth', 2);
            ylabel(ax1, 'amplitude', 'FontSize', opts.LabelFontSize);
            title(ax1, 'Envelope, IMFs, and inst. frequency', 'FontSize', opts.TitleFontSize);
            xlim(ax1, xl);
            set(ax1, 'XTickLabel', []);
            grid(ax1, 'on');

            ax2 = fig.Axes(2);
            EnvelopeMetricsPlotter.setAxisFontSize(ax2, opts.AxesFontSize);
            plot(ax2, t_env, imfs, 'LineWidth', 2);
            ylabel(ax2, 'amplitude', 'FontSize', opts.LabelFontSize);
            legendLabels = "imf" + (1:size(imfs, 1));
            legend(ax2, legendLabels, 'Location', 'northeast', 'FontSize', opts.LegendFontSize);
            xlim(ax2, xl);
            set(ax2, 'XTickLabel', []);
            grid(ax2, 'on');

            ax3 = fig.Axes(3);
            EnvelopeMetricsPlotter.setAxisFontSize(ax3, opts.AxesFontSize);
            plot(ax3, t_env, imfw, 'LineWidth', 2);
            xlabel(ax3, 'time (s)', 'FontSize', opts.LabelFontSize);
            ylabel(ax3, 'inst. freq. (Hz)', 'FontSize', opts.LabelFontSize);
            legendLabels = "imf" + (1:size(imfw, 1)) + " freq.";
            legend(ax3, legendLabels, 'Location', 'northeast', 'FontSize', opts.LegendFontSize);
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
            fig = EnvelopeMetricsPlotter.createSizedFigure([1 1], EnvelopeMetricsPlotter.ExteriorMargins, ...
                [], EnvelopeMetricsPlotter.FigureSizeInches);
        end

        function fig = createSizedFigure(panelSpec, extMargins, intMargins, sizeInches)
            % Passing Position at construction (rather than only setting it
            % afterward) skips stFig's WindowState='Maximized' step, which
            % does not reliably resize back down to an explicit Position in
            % a headless -batch session -- without this, the first figure
            % exported in a session could end up rendered at a different
            % physical size than later ones, throwing off the apparent
            % relative size of same-point-size fonts across figures.
            fig = stFig(panelSpec, extMargins, intMargins, 'Theme', 'light', 'Position', sizeInches);
            set(fig.Handle, 'Units', 'inches', 'Position', sizeInches);
        end

        function setAxisFontSize(ax, fontSize)
            % Tick label size; title/axis/legend labels set their own FontSize above.
            set(ax, 'FontSize', fontSize);
        end

    end

end
