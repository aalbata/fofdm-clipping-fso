function plot_block_diagram(figuresDir)
%PLOT_BLOCK_DIAGRAM Block diagram of the simulated baseband chain.
%   Solid blocks are used in every configuration. Dashed blocks are used only in the
%   configurations with subband filtering (filters) or clipping (limiter).
st = figure_style();
H = 1.50;                                          % figure height (in)
fig = new_figure(st, H);
ax = axes(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);
axis(ax, [0 st.width 0 H]);
axis(ax, 'off');
hold(ax, 'on');

bw = 0.50;  bh = 0.34;  gap = (st.width - 0.03 - 6 * bw) / 5;
x = 0.015 + (0:5) * (bw + gap);                    % left edges of the six columns
yTx = 1.08;  yRx = 0.08;
xc = x(6) + bw / 2;                                % center of the last column
yMul = 0.80;  yAdd = 0.60;  rc = 0.055;            % multiplier and adder

txLabels = {{'Random', 'bits'}, {'16-QAM', 'mapping'}, {'IFFT,', 'add CP'}, ...
    {'Subband', 'filter'}, {'Limiter'}, {'Scale \kappa'}};
txDashed = [false false false true true false];
for i = 1:6
    block(x(i), yTx, bw, bh, txLabels{i}, txDashed(i));
    if i < 6
        arrow([x(i) + bw, x(i + 1)], [yTx, yTx] + bh / 2);
    end
end
rxLabels = {{'Matched', 'filter'}, {'Remove', 'CP, FFT'}, {'Equalizer'}, {'16-QAM', 'demapping'}, {'Bit error', 'count'}};
rxDashed = [true false false false false];
for i = 1:5
    col = 7 - i;
    block(x(col), yRx, bw, bh, rxLabels{i}, rxDashed(i));
    if i < 5
        arrow([x(col), x(col - 1) + bw], [yRx, yRx] + bh / 2);
    end
end

% Channel: gain of the power record and additive noise.
arrow([xc, xc], [yTx, yMul + rc]);
node(xc, yMul, '\times');
arrow([xc, xc], [yMul - rc, yAdd + rc]);
node(xc, yAdd, '+');
arrow([xc, xc], [yAdd - rc, yRx + bh]);
block(x(3) + 0.10, yMul - 0.11, 2 * bw + gap - 0.10, 0.22, {'Measured power record'}, false);
arrow([x(4) + bw, xc - rc], [yMul, yMul]);
text(ax, (x(4) + bw + xc - rc) / 2, yMul + 0.06, 'a[m]', 'FontName', st.font, 'FontSize', st.fontSize, ...
    'FontAngle', 'italic', 'HorizontalAlignment', 'center');
arrow([x(5) + 0.30, xc - rc], [yAdd, yAdd]);
text(ax, x(5) + 0.27, yAdd, 'Noise', 'FontName', st.font, 'FontSize', st.fontSize, 'HorizontalAlignment', 'right');
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
