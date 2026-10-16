function fig = new_figure(st, height)
%NEW_FIGURE Invisible figure window of one column width and the given height in inches.
%   The size in inches is stored in the UserData property for SAVE_FIGURE.
fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 st.width height], 'Color', 'w', ...
    'UserData', [st.width height]);
end
