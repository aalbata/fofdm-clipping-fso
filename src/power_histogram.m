function [h, pSum, pCount, pMax] = power_histogram(p, xn, outRate, D, s2, edges)
%POWER_HISTOGRAM Histogram of the instantaneous power of a reconstructed transmitted waveform.
%   xn is a transmitted waveform zero-padded to a multiple of p.nPad samples (TRANSMIT_STREAM);
%   it is reconstructed at outRate times the sampling rate by band-limited interpolation. Only
%   the samples of the counted blocks, which start D + p.Nb sampling periods after the first
%   sample, are used. h counts the levels in decibels relative to s2 between the edges;
%   pSum, pCount, and pMax are the sum, number, and largest value of the powers.
xp = resample_stream(xn, outRate * p.nPad);
seg = xp(outRate * (D + p.Nb) + (1:outRate * p.M * p.Nb));
pw = abs(seg).^2;
lev = 10 * log10(pw / s2);
if any(lev >= edges(end))
    error('power_histogram:range', 'A peak exceeds the histogram range.');
end
h = histcounts(lev, edges);
pSum = sum(pw);
pCount = numel(pw);
pMax = max(pw);
end
