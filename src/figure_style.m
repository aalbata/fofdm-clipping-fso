function st = figure_style()
%FIGURE_STYLE Sizes, fonts, colors, markers, and labels shared by all figures.
%   The six configurations are, in this order: OFDM and F-OFDM without a limiter, and the
%   transmit chains clip, filter_clip, clip_filter, and rcf at their selected clipping ratios.
st.width = 3.54;                                   % figure width, one column (in)
st.wideWidth = 5.51;                               % figure width, one and a half columns (in)
st.font = 'Times New Roman';
st.fontSize = 7;                                   % pt
st.lineWidth = 0.7;                                % pt
st.axesLineWidth = 0.4;                            % pt
st.markerSize = 3;                                 % pt
st.names = {'ofdm', 'fofdm', 'clip', 'filter_clip', 'clip_filter', 'rcf'};
st.labels = {'OFDM', 'F-OFDM', 'Clipped OFDM', 'Filter-then-clip', 'Clip-then-filter', 'RCF'};
st.colors = [0.000 0.447 0.698; 0.835 0.369 0.000; 0.000 0.620 0.451; 0.800 0.475 0.655; ...
    0.902 0.624 0.000; 0.337 0.706 0.914];
st.markers = {'o', 's', '^', 'd', 'v', '>'};
st.reference = [0 0 0];                            % color of analytical reference curves
st.gray = [0.45 0.45 0.45];                        % color of additional curves
end
