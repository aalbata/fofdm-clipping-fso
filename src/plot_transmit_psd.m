function plot_transmit_psd(resultsDir, figuresDir)
%PLOT_TRANSMIT_PSD Power spectral density of the transmitted stream for five streams.
%   The density is symmetric about zero frequency; the upper part of the band and the
%   guard band are shown. OFDM and F-OFDM are drawn without limiter, F-OFDM also with the
%   clipping ratio p.gamma0 (gray), and clipped OFDM and joint clipping-assisted F-OFDM
%   with the selected clipping ratio.
p = sim_params();
st = figure_style();
T = readtable(fullfile(resultsDir, 'psd.csv'));
rows = T.freq_mhz >= 0;
f = T.freq_mhz(rows);

fig = new_figure(st, 1.6);
tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
hold(ax, 'on');
h = gobjects(1, 5);
h(1) = plot(ax, f, T.ofdm_db(rows), '-', 'Color', st.colors(1, :), 'LineWidth', st.lineWidth);
h(2) = plot(ax, f, T.fofdm_db(rows), '-', 'Color', st.colors(2, :), 'LineWidth', st.lineWidth);
h(3) = plot(ax, f, T.fofdm_gamma0_db(rows), '-', 'Color', st.gray, 'LineWidth', st.lineWidth);
h(4) = plot(ax, f, T.clipped_ofdm_db(rows), '-', 'Color', st.colors(3, :), 'LineWidth', st.lineWidth);
h(5) = plot(ax, f, T.joint_db(rows), '-', 'Color', st.colors(4, :), 'LineWidth', st.lineWidth);
set(ax, 'XLim', [10 p.Fs / 2e6], 'YLim', [-150 10], 'YTick', -150:25:0, 'XTick', 10:2:30);
style_axes(ax, st);
xlabel(ax, 'Frequency (MHz)');
ylabel(ax, 'Normalized PSD (dB)');
labels = {st.labels{1}, st.labels{2}, sprintf('F-OFDM, \\gamma_0 = %g', p.gamma0), st.labels{3}, st.labels{4}};
lg = legend(ax, h, labels, 'Location', 'southwest', 'FontSize', st.fontSize - 1, 'Box', 'off');
lg.ItemTokenSize = [14 18];
save_figure(fig, figuresDir, 'transmit_psd');
end
