function fig = new_figure(st, height, width)
%NEW_FIGURE Invisible figure window of the given height and width in inches.
%   The width defaults to one column. The size is stored in the UserData property for
%   SAVE_FIGURE.
if nargin < 3
    width = st.width;
end
fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 width height], 'Color', 'w', ...
    'UserData', [width height]);
end
