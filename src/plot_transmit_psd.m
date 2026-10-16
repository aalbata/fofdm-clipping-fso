function plot_transmit_psd(resultsDir, figuresDir)
%PLOT_TRANSMIT_PSD Power spectral density of the transmitted waveforms of the six configurations.
%   The density is symmetric about zero frequency; the upper part of the band and the
%   guard band are shown, normalized to the mean density for |f| <= p.inBandHz (psd.csv). The
%   shaded region, p.guardHz <= f <= Fs/2, is the upper part of the range over which the
%   guard-band level is averaged; RCF is dashed.
p = sim_params();
st = figure_style();
T = readtable(fullfile(resultsDir, 'psd.csv'));
rows = T.freq_mhz >= 0;
f = T.freq_mhz(rows);

fig = new_figure(st, 2.1);
tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
hold(ax, 'on');
patch(ax, [p.guardHz, p.Fs / 2, p.Fs / 2, p.guardHz] / 1e6, [-155 -155 10 10], 0.93 * [1 1 1], ...
    'EdgeColor', 'none');                          % range of the guard-band level
h = gobjects(1, numel(st.names));
for c = 1:numel(st.names)
    y = T.([st.names{c} '_db'])(rows);
    style = '-';
    if strcmp(st.names{c}, 'rcf')
        style = '--';                              % RCF coincides with OFDM in the guard band
    end
    h(c) = plot(ax, f, y, style, 'Color', st.colors(c, :), 'LineWidth', st.lineWidth);
end
set(ax, 'XLim', [10 p.Fs / 2e6], 'YLim', [-155 10], 'YTick', -150:25:0, 'XTick', 10:2:30, 'Layer', 'top');
style_axes(ax, st);
xlabel(ax, 'Frequency (MHz)');
ylabel(ax, 'Normalized PSD (dB)');
lg = legend(ax, h, st.labels, 'Location', 'southwest', 'FontSize', st.fontSize - 1, 'Box', 'off');
lg.ItemTokenSize = [14 18];
save_figure(fig, figuresDir, 'transmit_psd');
end
