function h = subband_filter(p)
%SUBBAND_FILTER Impulse response of the subband filter.
%   Windowed-sinc lowpass filter of p.L taps with unit gain at zero frequency. The passband
%   covers the occupied band (data subcarriers and the zero-frequency bin) plus
%   p.toneOffset subcarriers on each side.
n = (0:p.L-1).' - (p.L - 1) / 2;
W = (2 * max(p.k) + 1 + 2 * p.toneOffset) / p.N;   % passband width relative to the sampling rate
x = W * n;
pb = ones(size(x));
nz = x ~= 0;
pb(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
h = pb .* dolph_chebyshev(p.L, p.sidelobeDb);
h = h / sum(h);
end
