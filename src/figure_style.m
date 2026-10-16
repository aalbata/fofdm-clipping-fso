function st = figure_style()
%FIGURE_STYLE Sizes, fonts, colors, and markers shared by all figures.
%   The four configurations are, in this order: OFDM, F-OFDM, clipped OFDM, and joint
%   clipping-assisted F-OFDM.
st.width = 3.45;                                   % figure width, one column (in)
st.font = 'Times New Roman';
st.fontSize = 7;                                   % pt
st.lineWidth = 0.7;                                % pt
st.axesLineWidth = 0.4;                            % pt
st.markerSize = 3;                                 % pt
st.names = {'ofdm', 'fofdm', 'clipped_ofdm', 'joint'};
st.labels = {'OFDM', 'F-OFDM', 'Clipped OFDM', 'Joint clipping-assisted F-OFDM'};
st.colors = [0.000 0.447 0.698; 0.835 0.369 0.000; 0.000 0.620 0.451; 0.800 0.475 0.655];
st.markers = {'o', 's', '^', 'd'};
st.reference = [0 0 0];                            % color of analytical reference curves
st.gray = [0.45 0.45 0.45];                        % color of additional curves
end
