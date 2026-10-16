function plot_peak_ccdf(resultsDir, figuresDir)
%PLOT_PEAK_CCDF Exceedance probability of the instantaneous power of the transmitted waveforms.
%   Fraction of the samples of the waveform reconstructed at p.peakRate times the sampling
%   rate whose power exceeds the average power by more than the abscissa, for the six
%   configurations (ccdf.csv). The dashed line marks the exceedance p.peakProb that defines
%   the peak power.
p = sim_params();
st = figure_style();
cfg = configurations(p, resultsDir);
T = readtable(fullfile(resultsDir, 'ccdf.csv'), 'TextType', 'string');

fig = new_figure(st, 2.1);
tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
hold(ax, 'on');
h = gobjects(1, numel(cfg.name));
for c = 1:numel(cfg.name)
    rows = T.chain == cfg.chain{c} & (T.gamma == cfg.gamma(c) | (isinf(T.gamma) & isinf(cfg.gamma(c))));
    lev = T.level_db(rows);
    pr = T.exceedance(rows);
    h(c) = semilogy(ax, lev, pr, '-', 'Color', st.colors(c, :), 'LineWidth', st.lineWidth, ...
        'Marker', st.markers{c}, 'MarkerSize', st.markerSize, 'MarkerIndices', c:12:numel(lev), 'MarkerFaceColor', 'w');
end
yline(ax, p.peakProb, '--', 'Color', st.gray, 'LineWidth', st.lineWidth * 0.8);
set(ax, 'YScale', 'log', 'XLim', [2 12], 'YLim', [1e-6 1], 'YTick', 10.^(-6:0), 'XTick', 2:12);
style_axes(ax, st);
xlabel(ax, 'Instantaneous power relative to average power (dB)');
ylabel(ax, 'Exceedance probability');
lg = legend(ax, h, st.labels, 'Location', 'southwest', 'FontSize', st.fontSize - 1, 'Box', 'off');
lg.ItemTokenSize = [14 18];
save_figure(fig, figuresDir, 'peak_ccdf');
end
