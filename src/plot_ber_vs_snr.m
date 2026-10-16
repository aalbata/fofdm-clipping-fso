function plot_ber_vs_snr(resultsDir, figuresDir)
%PLOT_BER_VS_SNR BER of the four configurations under both transmitter models (sunny record).
%   Top panel: BER versus SNR, average-power model, with the closed form for 16-QAM
%   averaged over the record as a dotted line. Bottom panel: BER versus peak SNR,
%   peak-limited model; each curve is the simulated curve shifted by the peak-to-average
%   power ratio of the limiter output. The clipping ratio of each curve is taken from
%   summary.csv. A marker is drawn at every second simulated SNR value, at even values for
%   OFDM and clipped OFDM and at odd values for F-OFDM and joint clipping-assisted F-OFDM.
p = sim_params();
st = figure_style();
S = read_ber(resultsDir);
T = readtable(fullfile(resultsDir, 'summary.csv'), 'TextType', 'string');
V = readtable(fullfile(resultsDir, 'validation.csv'), 'TextType', 'string');
ref = clipping_reference(p.gammaGrid);
paprDb = 10 * log10(ref.papr);
models = ["avg", "peak"];
xLabels = {'SNR (dB)', 'PSNR (dB)'};
xLimits = {[4 24], [12 32]};
panels = {'(a)', '(b)'};
iSun = find(strcmp(p.records, 'sun'));

H = 2.95;                                          % figure height (in)
fig = new_figure(st, H);
% panels placed by hand, in inches from the lower left corner; the shared legend sits above panel (a)
left = 0.42;  bottom = [1.68 0.34];  height = 0.92;
for m = 1:2
    ax = axes(fig, 'Units', 'inches', 'Position', [left, bottom(m), st.width - left - 0.06, height]);
    hold(ax, 'on');
    h = gobjects(1, 4);
    for c = 1:4
        row = T.model == models(m) & T.configuration == st.names{c};
        ig = find(round(10 * p.gammaGrid) == round(10 * T.gamma(row)));
        ber = reshape(S.ber(iSun, T.filter(row) + 1, ig, :), 1, []);
        x = p.snrGrid;
        if m == 2
            x = x + paprDb(ig);
        end
        keep = ber > 0;
        h(c) = semilogy(ax, x(keep), ber(keep), '-', 'Color', st.colors(c, :), 'LineWidth', st.lineWidth, ...
            'Marker', st.markers{c}, 'MarkerSize', st.markerSize, 'MarkerIndices', (1 + mod(c + 1, 2)):2:nnz(keep), ...
            'MarkerFaceColor', 'w');
    end
    if m == 1
        rows = V.record == "sun";
        hRef = semilogy(ax, V.snr_db(rows), V.ber_theory(rows), ':', 'Color', st.reference, 'LineWidth', st.lineWidth);
        lg = legend(ax, [h, hRef], [st.labels, {'16-QAM closed form'}], 'Orientation', 'horizontal', ...
            'NumColumns', 2, 'FontSize', st.fontSize - 1, 'Box', 'off', 'Units', 'inches');
        lg.ItemTokenSize = [18 18];
    end
    set(ax, 'YScale', 'log', 'XLim', xLimits{m}, 'YLim', [1e-5 1], 'YTick', 10.^(-5:0));
    style_axes(ax, st);
    xlabel(ax, xLabels{m});
    ylabel(ax, 'BER');
    panel_label(ax, st, panels{m});
end
drawnow;
lg.Position(1:2) = [left + (st.width - left - 0.06 - lg.Position(3)) / 2, H - lg.Position(4)];
save_figure(fig, figuresDir, 'ber_vs_snr');
end
