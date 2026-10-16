function plot_block_diagram(figuresDir)
%PLOT_BLOCK_DIAGRAM Block diagram of the simulated baseband chain.
%   Solid blocks are used in every configuration. Dashed blocks are used only in the
%   transmit chains named below them: C (clip), FC (filter_clip), CF (clip_filter), and
%   RCF (rcf); the first subband filter is also used by F-OFDM, and the matched receive filter
%   with FC, CF, and F-OFDM.
st = figure_style();
W = st.wideWidth;
H = 1.95;                                          % figure height (in)
fig = new_figure(st, H, W);
ax = axes(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);
axis(ax, [0 W 0 H]);
axis(ax, 'off');
hold(ax, 'on');

nCol = 9;
bw = 0.52;  bh = 0.40;  gap = (W - 0.03 - nCol * bw) / (nCol - 1);
x = 0.015 + (0:nCol - 1) * (bw + gap);             % left edges of the columns
yTx = 1.40;  yRx = 0.22;
xc = x(nCol) + bw / 2;                             % center of the last column
yMul = 1.08;  yAdd = 0.86;  rc = 0.055;            % multiplier and adder

txLabels = {{'Random', 'bits'}, {'16-QAM', 'mapping'}, {'IFFT'}, {'Repeated', 'clipping,', 'filtering'}, ...
    {'Add CP'}, {'Subband', 'filter'}, {'Limiter', 'at 4{\itF}_s'}, {'Subband', 'filter'}, {'Scale \kappa'}};
txDashed = [false false false true false true true true false];
txTag = {'', '', '', 'RCF', '', 'FC, F-OFDM', 'C, FC, CF', 'CF', ''};
for i = 1:nCol
    block(x(i), yTx, bw, bh, txLabels{i}, txDashed(i));
    if ~isempty(txTag{i})
        text(ax, x(i) + bw / 2, yTx - 0.07, txTag{i}, 'FontName', st.font, 'FontSize', st.fontSize - 1, ...
            'FontAngle', 'italic', 'HorizontalAlignment', 'center');
    end
    if i < nCol
        arrow([x(i) + bw, x(i + 1)], [yTx, yTx] + bh / 2);
    end
end
rxLabels = {{'Matched', 'filter'}, {'Remove', 'CP, FFT'}, {'Equalizer'}, {'16-QAM', 'demapping'}, {'Bit error', 'count'}};
rxDashed = [true false false false false];
for i = 1:5
    col = nCol + 1 - i;
    block(x(col), yRx, bw, bh, rxLabels{i}, rxDashed(i));
    if i < 5
        arrow([x(col), x(col - 1) + bw], [yRx, yRx] + bh / 2);
    end
end
text(ax, x(nCol) + bw, yRx - 0.08, 'FC, CF, F-OFDM', 'FontName', st.font, 'FontSize', st.fontSize - 1, ...
    'FontAngle', 'italic', 'HorizontalAlignment', 'right');

% Channel: amplitude gain of the power record and additive noise.
arrow([xc, xc], [yTx, yMul + rc]);
node(xc, yMul, '\times');
arrow([xc, xc], [yMul - rc, yAdd + rc]);
node(xc, yAdd, '+');
arrow([xc, xc], [yAdd - rc, yRx + bh]);
block(x(5), yMul - 0.11, 2 * bw + gap, 0.22, {'Measured power record'}, false);
arrow([x(6) + bw, xc - rc], [yMul, yMul]);
text(ax, (x(6) + bw + xc - rc) / 2, yMul + 0.07, '{\ita}[{\itm}]', 'FontName', st.font, 'FontSize', st.fontSize, ...
    'HorizontalAlignment', 'center');
arrow([x(7) + 0.40, xc - rc], [yAdd, yAdd]);
text(ax, x(7) + 0.37, yAdd, 'Noise', 'FontName', st.font, 'FontSize', st.fontSize, 'HorizontalAlignment', 'right');
text(ax, x(1), yTx - 0.10, 'Transmitter', 'FontName', st.font, 'FontSize', st.fontSize, 'FontAngle', 'italic');
text(ax, x(1), yRx + bh / 2, 'Receiver', 'FontName', st.font, 'FontSize', st.fontSize, 'FontAngle', 'italic');

save_figure(fig, figuresDir, 'block_diagram');

    function block(x0, y0, w, h, label, dashed)
        %BLOCK Rectangle with a centered label; dashed outline for optional blocks.
        style = '-';
        if dashed
            style = '--';
        end
        rectangle(ax, 'Position', [x0 y0 w h], 'LineWidth', 0.5, 'LineStyle', style, 'EdgeColor', 'k');
        text(ax, x0 + w / 2, y0 + h / 2, label, 'FontName', st.font, 'FontSize', st.fontSize, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    function node(x0, y0, symbol)
        %NODE Circle with an operator symbol.
        rectangle(ax, 'Position', [x0 - rc, y0 - rc, 2 * rc, 2 * rc], 'Curvature', [1 1], 'LineWidth', 0.5);
        text(ax, x0, y0, symbol, 'FontName', st.font, 'FontSize', st.fontSize, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    function arrow(xs, ys)
        %ARROW Straight line from the first to the second point with an arrowhead at the end.
        d = [xs(2) - xs(1), ys(2) - ys(1)];
        d = d / norm(d);
        n = [-d(2), d(1)];
        tip = [xs(2), ys(2)];
        len = 0.045;  half = 0.018;
        base = tip - len * d;
        plot(ax, [xs(1), base(1)], [ys(1), base(2)], 'k-', 'LineWidth', 0.5);
        patch(ax, [tip(1), base(1) + half * n(1), base(1) - half * n(1)], ...
            [tip(2), base(2) + half * n(2), base(2) - half * n(2)], 'k', 'EdgeColor', 'none');
    end
end
