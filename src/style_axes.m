function style_axes(ax, st)
%STYLE_AXES Apply the shared font, line widths, grid, and box to an axes.
%   Negative tick labels are written with a minus sign instead of a hyphen.
set(ax, 'FontName', st.font, 'FontSize', st.fontSize, 'LineWidth', st.axesLineWidth, 'Box', 'on', ...
    'TickDir', 'in', 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'GridAlpha', 0.15, 'XColor', 'k', 'YColor', 'k');
ax.XTick = ax.XTick;                               % positions fixed together with their labels
ax.YTick = ax.YTick;
ax.XTickLabel = strrep(ax.XTickLabel, '-', char(8722));
ax.YTickLabel = strrep(ax.YTickLabel, '-', char(8722));
end
