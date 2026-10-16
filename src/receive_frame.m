function R = receive_frame(p, z, Dc)
%RECEIVE_FRAME FFT outputs on the data subcarriers for the counted blocks of a frame.
%   z is the stream at the FFT input and Dc the delay of the filter cascade in samples.
%   The FFT window of each block starts p.delta samples before the end of its cyclic prefix.
%   R is Nd-by-M with the unitary transform scaling.
n0 = (1:p.M) * p.Nb + p.Ncp - p.delta + Dc;        % samples before the first sample of each window
idx = n0 + (1:p.N).';                              % N-by-M sample indices (1-based)
R = fft(z(idx), [], 1) / sqrt(p.N);
R = R(p.bins, :);
end
