function plot_tradeoff(resultsDir, figuresDir)
%PLOT_TRADEOFF Guard-band level versus required peak SNR over the clipping ratio.
%   Sunny record, BER target p.selectionTarget, peak-limited model. Each line connects the
%   clipping ratios of one transmit chain from its selected clipping ratio (filled marker) to
%   3.0 in steps of 0.1. OFDM and F-OFDM without a limiter are shown as single markers.
p = sim_params();
st = figure_style();
cfg = configurations(p, resultsDir);
T = readtable(fullfile(resultsDir, 'required_snr.csv'), 'TextType', 'string');
[~, ~, cal] = read_calibration(p, resultsDir);
target = p.selectionTarget;

fig = new_figure(st, 2.4);
tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
hold(ax, 'on');
h = gobjects(1, numel(cfg.name));
for c = 1:numel(cfg.name)
    color = st.colors(c, :);
    if isinf(cfg.gamma(c))
        [x, y] = point(cfg.chain{c}, Inf);
        h(c) = plot(ax, x, y, st.markers{c}, 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', st.markerSize + 1);
        continue
    end
    g = p.gammaGrid(p.gammaGrid >= cfg.gamma(c) - 1e-9 & ~isinf(p.gammaGrid));
    x = zeros(size(g));
    y = zeros(size(g));
    for i = 1:numel(g)
        [x(i), y(i)] = point(cfg.chain{c}, g(i));
    end
    h(c) = plot(ax, x, y, '-', 'Color', color, 'LineWidth', st.lineWidth, 'Marker', st.markers{c}, ...
        'MarkerSize', st.markerSize - 1);
end
% selected clipping ratios, larger markers first so that coinciding points stay visible
big = [0 0 st.markerSize + 1 st.markerSize + 3 st.markerSize + 1 st.markerSize + 1];
for c = [4 3 5 6]                                % filter-then-clip before clipped OFDM
    [xs, ys] = point(cfg.chain{c}, cfg.gamma(c));
    plot(ax, xs, ys, st.markers{c}, 'Color', st.colors(c, :), 'MarkerFaceColor', st.colors(c, :), 'MarkerSize', big(c));
end
set(ax, 'XLim', [20 26], 'YLim', [-140 -10], 'XTick', 20:26, 'YTick', -140:20:-20);
style_axes(ax, st);
xlabel(ax, 'Required PSNR (dB)');
ylabel(ax, 'Guard-band level (dB)');
lg = legend(ax, h, st.labels, 'Location', 'west', 'FontSize', st.fontSize - 1, 'Box', 'off');
lg.ItemTokenSize = [14 18];
save_figure(fig, figuresDir, 'tradeoff');

    function [x, y] = point(chain, gamma)
        %POINT Required peak SNR and guard-band level of a chain at a clipping ratio.
        rr = T.record == "sun" & T.chain == chain & T.ber_target == target & ...
            (T.gamma == gamma | (isinf(T.gamma) & isinf(gamma)));
        rc = cal.chain == chain & (cal.gamma == gamma | (isinf(cal.gamma) & isinf(gamma)));
        x = T.psnr_req_db(rr);
        y = cal.guard_band_db(rc);
    end
end
