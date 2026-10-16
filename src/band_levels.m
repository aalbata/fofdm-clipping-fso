function [psd, freq, guardDb, outsideDb] = band_levels(p, acc, rate)
%BAND_LEVELS Normalized spectral density within half the sampling rate and its band levels.
%   acc is the sum of periodograms of segments of p.psdSegment sampling periods taken at
%   rate times the sampling rate. psd is the density at the p.psdSegment frequencies freq
%   (Hz) within half the sampling rate, relative to the mean density for |f| <= p.inBandHz.
%   guardDb is the mean of psd over p.guardHz <= |f| <= Fs/2 (both band
%   edges), in decibels. outsideDb is the power beyond half the sampling rate relative to
%   the total power in decibels (NaN for rate = 1).
n = p.psdSegment;
full = fftshift(acc);
center = numel(full) / 2 + (-n / 2 + 1:n / 2).';
freq = (-n / 2:n / 2 - 1).' * p.Fs / n;
psd = full(center) / mean(full(center(abs(freq) <= p.inBandHz)));
guardDb = 10 * log10(mean(psd(abs(freq) >= p.guardHz & freq > -p.Fs / 2)));   % same number of bins at both edges
outsideDb = NaN;
if rate > 1
    outsideDb = 10 * log10(1 - sum(full(center)) / sum(full));
end
end
