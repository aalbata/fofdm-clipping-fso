function y = filter_stream(x, h)
%FILTER_STREAM Full convolution of a column stream with a filter, computed with the FFT.
%   The output has numel(x) + numel(h) - 1 samples.
n = numel(x) + numel(h) - 1;
nfft = 2^nextpow2(n);
y = ifft(fft(x, nfft) .* fft(h, nfft));
y = y(1:n);
end
