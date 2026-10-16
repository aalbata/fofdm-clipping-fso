function s = required_snr(snr, ber, target, method)
%REQUIRED_SNR SNR at which a BER curve crosses a target.
%   snr is an increasing grid in dB and ber the BER at each grid point. The crossing is
%   located between the first pair of adjacent grid points i and i+1 with
%   ber(i) >= target > ber(i+1) > 0, on the parabola through log10(ber) at three grid
%   points: i-1, i, i+1 when point i-1 exists and has nonzero BER, otherwise i, i+1, i+2 when
%   point i+2 exists and has nonzero BER. Without a third point, or when the parabola has no crossing
%   between the two grid points, the crossing is taken on the straight line through
%   log10(ber) at i and i+1. The result is NaN when no such pair exists. With method
%   'linear' the crossing is always taken on the straight line; the default is 'parabola'.
if nargin < 4
    method = 'parabola';
end
s = NaN;
n = numel(snr);
for i = 1:n - 1
    if ber(i) >= target && ber(i + 1) < target && ber(i + 1) > 0
        y = log10(target);
        frac = (log10(ber(i)) - y) / (log10(ber(i)) - log10(ber(i + 1)));
        s = snr(i) + frac * (snr(i + 1) - snr(i));
        if strcmp(method, 'linear')
            return
        end
        if i > 1 && ber(i - 1) > 0
            j = i - 1:i + 1;
        elseif i + 2 <= n && ber(i + 2) > 0
            j = i:i + 2;
        else
            return
        end
        c = polyfit(snr(j) - snr(i), log10(ber(j)), 2);
        c(3) = c(3) - y;
        root = roots(c);
        root = real(root(abs(imag(root)) < 1e-12)) + snr(i);
        root = root(root >= snr(i) & root <= snr(i + 1));
        if isscalar(root)
            s = root;
        end
        return
    end
end
end
