function plot_clipping_ratio_sweep(resultsDir, figuresDir)
%PLOT_CLIPPING_RATIO_SWEEP BER and required peak SNR versus the clipping ratio (sunny record).
%   Top panel: BER of clipped OFDM (lines) and joint clipping-assisted F-OFDM (markers)
%   versus the clipping ratio at three SNR values, average-power model.
%   Bottom panel: required peak SNR versus the clipping ratio for the BER targets,
%   peak-limited model, with the minimum of each clipped-OFDM curve marked by a filled
%   circle. Clipped OFDM is drawn with lines and joint clipping-assisted F-OFDM with markers.
%   Values above the upper limit of the vertical axis are not drawn.
p = sim_params();
st = figure_style();
S = read_ber(resultsDir);
[req, ~, paprDb] = required_table(p, S);
finite = ~isinf(p.gammaGrid);
gamma = p.gammaGrid(finite);
snrShown = [12 16 20];
shade = [0.55 0.30 0.00];                          % gray level of the three SNR curves
iSun = find(strcmp(p.records, 'sun'));

fig = new_figure(st, 2.45);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
hold(ax, 'on');
h = gobjects(1, numel(snrShown));
for i = 1:numel(snrShown)
    iSnr = find(p.snrGrid == snrShown(i));
    color = shade(i) * [1 1 1];
    h(i) = semilogy(ax, gamma, squeeze(S.ber(iSun, 1, finite, iSnr)), '-', 'Color', color, 'LineWidth', st.lineWidth);
    semilogy(ax, gamma, squeeze(S.ber(iSun, 2, finite, iSnr)), 'LineStyle', 'none', 'Marker', st.markers{4}, ...
        'Color', color, 'MarkerSize', st.markerSize - 0.5, 'LineWidth', 0.4);
end
set(ax, 'YScale', 'log', 'XLim', [0.2 3], 'YLim', [1e-5 1], 'YTick', 10.^(-5:0));
style_axes(ax, st);
ylabel(ax, 'BER');
lg = legend(ax, h, arrayfun(@(s) sprintf('SNR = %d dB', s), snrShown, 'UniformOutput', false), ...
    'Location', 'southwest', 'FontSize', st.fontSize - 1, 'Box', 'off');
lg.ItemTokenSize = [12 14];
panel_label(ax, st, '(a)');

ax = nexttile(tl);
hold(ax, 'on');
lineStyles = {'-', '--'};
h = gobjects(1, numel(p.berTargets));
labels = cell(1, numel(p.berTargets));
for it = 1:numel(p.berTargets)
    psnr = reshape(req(iSun, :, finite, it), 2, []) + paprDb(finite);
    h(it) = plot(ax, gamma, psnr(1, :), lineStyles{it}, 'Color', st.colors(3, :), 'LineWidth', st.lineWidth);
    plot(ax, gamma, psnr(2, :), 'LineStyle', 'none', 'Marker', st.markers{4}, 'Color', st.colors(4, :), ...
        'MarkerSize', st.markerSize - 0.5, 'LineWidth', 0.4);
    [best, iBest] = min(psnr(1, :));
    plot(ax, gamma(iBest), best, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', st.markerSize - 0.5);
    labels{it} = sprintf('BER = 10^{%s%d}', char(8722), abs(round(log10(p.berTargets(it)))));
end
set(ax, 'XLim', [0.2 3], 'YLim', [16 30]);
style_axes(ax, st);
xlabel(ax, 'Clipping ratio \gamma');
ylabel(ax, 'Required PSNR (dB)');
lg = legend(ax, h, labels, 'Location', 'southwest', 'FontSize', st.fontSize, 'Box', 'off');
lg.ItemTokenSize = [14 18];
panel_label(ax, st, '(b)');

save_figure(fig, figuresDir, 'clipping_ratio_sweep');
end
