function s = rcf_stream(p, X, A)
%RCF_STREAM Block stream after repeated clipping and frequency-domain filtering.
%   X is the Nd-by-B array of data symbols and A the clipping threshold. Each OFDM symbol
%   is interpolated by the factor p.os (zero padding of the spectrum), passed through the
%   envelope limiter, and transformed back; all spectral components outside the data
%   subcarriers are then set to zero. The clipping and filtering are repeated
%   p.rcfIterations times with the same threshold. The cyclic prefix is added to the final
%   symbols, which are not clipped again. s is the column stream of B blocks at the
%   sampling rate, with the scaling of OFDM_FRAME.
B = size(X, 2);
n = p.os * p.N;
bins = mod(p.k, n) + 1;                            % data subcarriers on the interpolated grid
Y = X;
for it = 1:p.rcfIterations
    F = zeros(n, B);
    F(bins, :) = Y;
    x = ifft(F, [], 1) * (n / sqrt(p.N));          % symbol interpolated by p.os, scaling of (1)
    x = envelope_limiter(x, A);
    F = fft(x, [], 1) * (sqrt(p.N) / n);
    Y = F(bins, :);                                % components outside the data subcarriers removed
end
Xf = zeros(p.N, B);
Xf(p.bins, :) = Y;
x = ifft(Xf, [], 1) * sqrt(p.N);
blocks = [x(p.N - p.Ncp + 1:p.N, :); x];
s = blocks(:);
end
