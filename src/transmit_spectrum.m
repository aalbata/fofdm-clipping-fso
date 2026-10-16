function transmit_spectrum(resultsDir)
%TRANSMIT_SPECTRUM Spectral density, guard-band level, power beyond half the sampling rate, and peak
%   excess of the transmitted stream.
%   Streams, for no filter and for the subband filter:
%   ofdm, fofdm                         no limiter
%   ofdm_gamma0, fofdm_gamma0           clipping ratio p.gamma0, limiter at p.os times the
%                                       sampling rate
%   clipped_ofdm, joint                 clipping ratio selected in selected_gamma.csv (sunny
%                                       record), limiter at p.os times the sampling rate
%   clipped_ofdm_rate1, joint_rate1     selected clipping ratio, limiter at the sampling rate
%   clipped_ofdm_lowpass, joint_lowpass selected clipping ratio, limiter at p.os times the
%                                       sampling rate, output restricted to the sampling
%                                       bandwidth by an ideal lowpass filter
%   Writes to resultsDir:
%   psd.csv             density of the first six streams for frequencies within half the
%                       sampling rate, in decibels relative to the mean in-band density
%                       (mean over frequencies up to 20 MHz)
%   transmit_stats.csv  for every stream: mean density over the outer guard band (24 MHz to
%                       half the sampling rate) relative to the mean in-band density, power
%                       beyond half the sampling rate
%                       relative to the total power, and the median and maximum over the
%                       blocks of the peak magnitude of the stream interpolated to 16
%                       times the sampling rate, relative to the clipping threshold
%   Estimator: averaged periodograms of Kaiser-windowed segments (shape parameter 30) of
%   8192 sampling periods with a hop of 2048 sampling periods, taken from each frame
%   without its first and last block. The density uses 60 frames and the peaks the first
%   6 of them.
p = sim_params();
ctx = sim_context(p);
sel = readtable(fullfile(resultsDir, 'selected_gamma.csv'), 'TextType', 'string');
gammaSel = [sel.gamma(sel.record == "sun" & sel.filter == 0), sel.gamma(sel.record == "sun" & sel.filter == 1)];
names = {'ofdm', 'fofdm', 'ofdm_gamma0', 'fofdm_gamma0', 'clipped_ofdm', 'joint', ...
    'clipped_ofdm_rate1', 'joint_rate1', 'clipped_ofdm_lowpass', 'joint_lowpass'};
filt = [0 1 0 1 0 1 0 1 0 1];
gamma = [Inf, Inf, p.gamma0, p.gamma0, gammaSel, gammaSel, gammaSel];
limiterRate = [1 1 p.os p.os p.os p.os 1 1 p.os p.os];   % sampling-rate multiple at the limiter
streamRate = [1 1 p.os p.os p.os p.os 1 1 1 1];          % sampling-rate multiple of the analyzed stream
nC = numel(names);
nPlot = 6;                                         % streams written to psd.csv

nSeg = 8192;                                       % segment length in sampling periods
nFrames = 60;
nPeakFrames = 6;
peakRate = 16;                                     % sampling-rate multiple used to measure peaks
D = [0, p.L / 2];                                  % transmit delay without and with the subband filter
acc = cell(1, nC);
count = zeros(1, nC);
peaks = cell(1, nC);
win = cell(1, nC);
for c = 1:nC
    acc{c} = zeros(streamRate(c) * nSeg, 1);
    win{c} = kaiser_window(streamRate(c) * nSeg, 30);
end
stream = RandStream('mt19937ar', 'Seed', 1001);
for f = 1:nFrames
    bits = randi(stream, [0 1], 4 * p.Nd, p.B);
    s0 = ofdm_frame(p, bits);
    u = {s0, filter_stream(s0, ctx.h)};
    for c = 1:nC
        T = numel(u{filt(c) + 1});
        x = [u{filt(c) + 1}; zeros(p.nPad - T, 1)];
        if ~isinf(gamma(c))
            A = gamma(c) * sqrt(ctx.sigma2(filt(c) + 1));
            x = envelope_limiter(resample_stream(x, limiterRate(c) * p.nPad), A);
            x = resample_stream(x, streamRate(c) * p.nPad);
            if f <= nPeakFrames
                xp = abs(resample_stream(x, peakRate * p.nPad));
                first = peakRate * (D(filt(c) + 1) + (1:p.M) * p.Nb);
                blockPeak = max(xp(first + (1:peakRate * p.Nb).'), [], 1);
                peaks{c} = [peaks{c}, 20 * log10(blockPeak / A)];
            end
        end
        n = streamRate(c) * nSeg;
        t = x(streamRate(c) * p.Nb + 1:streamRate(c) * (T - p.Nb));
        for start = 0:n / 4:numel(t) - n
            acc{c} = acc{c} + abs(fft(win{c} .* t(start + (1:n)))).^2;
            count(c) = count(c) + 1;
        end
    end
end

freq = (-nSeg / 2:nSeg / 2 - 1).' * p.Fs / nSeg;   % frequencies within half the sampling rate
inBand = abs(freq) <= 20e6;
guard = abs(freq) >= 24e6;
psd = zeros(nSeg, nC);
guardDb = zeros(1, nC);
outsideDb = NaN(1, nC);
for c = 1:nC
    full = fftshift(acc{c} / count(c));
    center = numel(full) / 2 + (-nSeg / 2 + 1:nSeg / 2).';
    psd(:, c) = full(center) / mean(full(center(inBand)));
    guardDb(c) = 10 * log10(mean(psd(guard, c)));
    if streamRate(c) > 1
        outsideDb(c) = 10 * log10(1 - sum(full(center)) / sum(full));
    end
end

fid = fopen(fullfile(resultsDir, 'psd.csv'), 'w');
fprintf(fid, 'freq_mhz');
fprintf(fid, ',%s_db', names{1:nPlot});
fprintf(fid, '\n');
fprintf(fid, ['%.4f' repmat(',%.3f', 1, nPlot) '\n'], [freq / 1e6, round(1e4 * log10(psd(:, 1:nPlot))) / 1e3 + 0].');
fclose(fid);

fid = fopen(fullfile(resultsDir, 'transmit_stats.csv'), 'w');
fprintf(fid, ['configuration,filter,gamma,limiter_rate,stream_rate,guard_band_db,outside_band_db,' ...
    'peak_excess_median_db,peak_excess_max_db\n']);
for c = 1:nC
    if isinf(gamma(c))
        fprintf(fid, '%s,%d,%s,NaN,%d,%s,NaN,NaN,NaN\n', names{c}, filt(c), gamma_text(gamma(c)), ...
            streamRate(c), num_text(guardDb(c), 2));
    else
        fprintf(fid, '%s,%d,%s,%d,%d,%s,%s,%s,%s\n', names{c}, filt(c), gamma_text(gamma(c)), ...
            limiterRate(c), streamRate(c), num_text(guardDb(c), 2), num_text(outsideDb(c), 2), ...
            num_text(median(peaks{c}), 2), num_text(max(peaks{c}), 2));
    end
end
fclose(fid);
end

function w = kaiser_window(n, beta)
%KAISER_WINDOW Symmetric Kaiser window of n samples with shape parameter beta.
x = (2 * (0:n - 1).' - (n - 1)) / (n - 1);
w = besseli(0, beta * sqrt(1 - x.^2)) / besseli(0, beta);
end
