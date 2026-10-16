function [acc, count] = welch_sum(acc, count, t, win)
%WELCH_SUM Add the periodograms of the windowed segments of a stream to a running sum.
%   Segments have the length of the window win and a hop of a quarter of that length.
%   acc is the running sum of squared FFT magnitudes and count the number of segments.
n = numel(win);
for start = 0:n / 4:numel(t) - n
    acc = acc + abs(fft(win .* t(start + (1:n)))).^2;
    count = count + 1;
end
end
