function y = resample_stream(x, nOut)
%RESAMPLE_STREAM Band-limited change of the sampling rate of a column stream by the FFT.
%   y has nOut samples and represents the same time span as x, taken as one period of a
%   periodic signal. For nOut > numel(x) the stream is interpolated; for nOut < numel(x)
%   it is restricted to the lower bandwidth (ideal lowpass filter) and decimated. numel(x)
%   and nOut are even. The component at half the lower sampling rate is shared equally
%   between the positive and the negative frequency, so that interpolation followed by
%   decimation returns the original stream.
nIn = numel(x);
if nOut == nIn
    y = x;
    return
end
X = fft(x);
h = min(nIn, nOut) / 2;                            % index of half the lower sampling rate
Y = zeros(nOut, 1);
Y(1:h) = X(1:h);
Y(nOut - h + 2:nOut) = X(nIn - h + 2:nIn);
if nOut > nIn
    Y(h + 1) = X(h + 1) / 2;
    Y(nOut - h + 1) = X(h + 1) / 2;
else
    Y(h + 1) = X(h + 1) + X(nIn - h + 1);
end
y = ifft(Y) * (nOut / nIn);
end
