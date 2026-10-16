function save_figure(fig, figuresDir, name)
%SAVE_FIGURE Write a figure as vector PDF and as PNG at 600 dots per inch, then close it.
%   fig is a figure created by NEW_FIGURE. Both files have the size stored by NEW_FIGURE.
sizeIn = fig.UserData;
exportgraphics(fig, fullfile(figuresDir, [name '.pdf']), 'ContentType', 'vector', 'Padding', 'figure', ...
    'Width', sizeIn(1), 'Height', sizeIn(2), 'Units', 'inches');
exportgraphics(fig, fullfile(figuresDir, [name '.png']), 'Resolution', 600, 'Padding', 'figure', ...
    'Width', sizeIn(1), 'Height', sizeIn(2), 'Units', 'inches');
close(fig);
end
