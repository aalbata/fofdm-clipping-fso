function plot_ber_vs_snr(resultsDir, figuresDir)
%PLOT_BER_VS_SNR BER of the six configurations under both transmitter models (sunny record).
%   Top panel: BER versus SNR, average-power model, with the closed form for 16-QAM averaged
%   over the record as a dotted line. Bottom panel: BER versus peak SNR, peak-limited model;
%   each curve is the simulated curve shifted by the PAPR of its transmitted waveform.
%   Points with no bit errors are not drawn. Markers are drawn every 2 dB of SNR, offset by
%   0.25 dB between consecutive configurations.
p = sim_params();
st = figure_style();
S = read_ber(resultsDir);
cfg = configurations(p, resultsDir);
[~, ~, cal] = read_calibration(p, resultsDir);
paprDb = chain_papr(p, cal);
V = readtable(fullfile(resultsDir, 'validation.csv'), 'TextType', 'string');
xLabels = {'SNR (dB)', 'PSNR (dB)'};
xLimits = {[4 22], [12 30]};
panels = {'(a)', '(b)'};
iSun = find(strcmp(p.records, 'sun'));
nCfg = numel(cfg.name);

H = 3.2;                                           % figure height (in)
fig = new_figure(st, H);
left = 0.42;  bottom = [1.68 0.34];  height = 0.96;
for m = 1:2
    ax = axes(fig, 'Units', 'inches', 'Position', [left, bottom(m), st.width - left - 0.06, height]);
    hold(ax, 'on');
    h = gobjects(1, nCfg);
    for c = 1:nCfg
        ber = reshape(S.ber(iSun, cfg.ic(c), cfg.ig(c), :), 1, []);
        x = p.snrGrid;
        if m == 2
            x = x + paprDb(cfg.ic(c), cfg.ig(c));
        end
        keep = find(ber > 0);
        marks = find(abs(mod(p.snrGrid(keep) - 0.25 * (c - 1), 2)) < 1e-9);   % every 2 dB, offset by 0.25 dB
        h(c) = semilogy(ax, x(keep), ber(keep), '-', 'Color', st.colors(c, :), 'LineWidth', st.lineWidth, ...
            'Marker', st.markers{c}, 'MarkerSize', st.markerSize, 'MarkerIndices', marks, ...
            'MarkerFaceColor', 'w');
    end
    if m == 1
        rows = V.record == "sun";
        hRef = semilogy(ax, V.snr_db(rows), V.ber_theory(rows), ':', 'Color', st.reference, 'LineWidth', st.lineWidth);
        lg = legend(ax, [h, hRef], [st.labels, {'16-QAM closed form'}], 'Orientation', 'horizontal', ...
            'NumColumns', 3, 'FontSize', st.fontSize - 1, 'Box', 'off', 'Units', 'inches');
        lg.ItemTokenSize = [18 18];
    end
    set(ax, 'YScale', 'log', 'XLim', xLimits{m}, 'YLim', [1e-5 1], 'YTick', 10.^(-5:0), 'XTick', xLimits{m}(1):2:xLimits{m}(2));
    style_axes(ax, st);
    xlabel(ax, xLabels{m});
    ylabel(ax, 'BER');
    panel_label(ax, st, panels{m});
end
drawnow;
lg.Position(1:2) = [left + (st.width - left - 0.06 - lg.Position(3)) / 2, H - lg.Position(4)];
save_figure(fig, figuresDir, 'ber_vs_snr');
end
