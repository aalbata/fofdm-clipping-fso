function plot_clipping_ratio_sweep(resultsDir, figuresDir)
%PLOT_CLIPPING_RATIO_SWEEP Required SNR, PAPR, and required peak SNR versus the clipping ratio.
%   Sunny record, BER target p.selectionTarget. (a) Required SNR under the average-power
%   model. (b) PAPR of the transmitted waveform (peak power at an exceedance of p.peakProb,
%   waveform reconstructed at p.peakRate times the sampling rate), with the PAPR (8) of the
%   limiter samples as a black dashed line. (c) Required peak SNR under the peak-limited
%   model, with the selected clipping ratio of each chain as a large filled marker; the marker
%   of filter_clip is drawn first and larger, so that the coinciding marker of clip stays
%   visible. Horizontal dotted lines show OFDM and F-OFDM without a limiter.
p = sim_params();
st = figure_style();
T = readtable(fullfile(resultsDir, 'required_snr.csv'), 'TextType', 'string');
sel = readtable(fullfile(resultsDir, 'selected_gamma.csv'), 'TextType', 'string');
ref = clipping_reference(p.gammaGrid(1:end - 1));
target = p.selectionTarget;
gamma = p.gammaGrid(1:end - 1);
cols = {'snr_req_db', 'papr_db', 'psnr_req_db'};
yLabels = {'Required SNR (dB)', 'PAPR (dB)', 'Required PSNR (dB)'};
yLimits = {[15 23], [0 12], [18 28]};
panels = {'(a)', '(b)', '(c)'};
unclipped = {{'clip', 1}, {'filter_clip', 2}};       % chain and color index of OFDM and F-OFDM

fig = new_figure(st, 3.9);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:3
    ax = nexttile(tl);
    hold(ax, 'on');
    h = gobjects(1, numel(p.chains));
    ySel = zeros(1, numel(p.chains));
    gSel = zeros(1, numel(p.chains));
    for ic = 1:numel(p.chains)
        rows = T.record == "sun" & T.chain == p.chains{ic} & T.ber_target == target & ~isinf(T.gamma);
        [g, order] = sort(T.gamma(rows));
        y = T.(cols{m})(rows);
        y = y(order);
        color = st.colors(2 + ic, :);
        h(ic) = plot(ax, g, y, '-', 'Color', color, 'LineWidth', st.lineWidth, 'Marker', st.markers{2 + ic}, ...
            'MarkerSize', st.markerSize - 1, 'MarkerIndices', 1 + mod(ic, 2):2:numel(g));
        gSel(ic) = sel.gamma(sel.chain == p.chains{ic});
        ySel(ic) = y(abs(g - gSel(ic)) < 1e-9);
    end
    if m == 3
        % selected clipping ratios, larger markers first so that coinciding points stay visible
        big = st.markerSize + [1 3 1 1];
        [~, drawOrder] = sort(big, 'descend');
        for ic = drawOrder
            color = st.colors(2 + ic, :);
            plot(ax, gSel(ic), ySel(ic), st.markers{2 + ic}, 'Color', color, 'MarkerFaceColor', color, ...
                'MarkerSize', big(ic));
        end
    end
    for u = 1:2
        row = T.record == "sun" & T.chain == unclipped{u}{1} & T.ber_target == target & isinf(T.gamma);
        yline(ax, T.(cols{m})(row), ':', 'Color', st.colors(unclipped{u}{2}, :), 'LineWidth', st.lineWidth);
    end
    if m == 2
        hRef = plot(ax, gamma, 10 * log10(ref.papr), '--', 'Color', st.reference, 'LineWidth', st.lineWidth * 0.8);
    end
    set(ax, 'XLim', [1 3], 'YLim', yLimits{m}, 'XTick', 1:0.2:3);
    style_axes(ax, st);
    ylabel(ax, yLabels{m});
    if m == 3
        xlabel(ax, 'Clipping ratio \gamma');
    end
    panel_label(ax, st, panels{m});
    if m == 1
        lg = legend(ax, h, st.labels(3:6), 'Location', 'northeast', 'NumColumns', 1, ...
            'FontSize', st.fontSize - 1, 'Box', 'off');
        lg.ItemTokenSize = [14 18];
    elseif m == 2
        lg = legend(ax, hRef, {'Limiter samples, (8)'}, 'Location', 'southeast', 'FontSize', st.fontSize - 1, 'Box', 'off');
        lg.ItemTokenSize = [14 18];
    end
end
save_figure(fig, figuresDir, 'clipping_ratio_sweep');
end
