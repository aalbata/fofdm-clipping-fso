function panel_label(ax, st, label)
%PANEL_LABEL Write a panel label such as (a) in the upper right corner of an axes.
text(ax, 0.97, 0.90, label, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontName', st.font, 'FontSize', st.fontSize);
end
