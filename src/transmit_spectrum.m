function transmit_spectrum(resultsDir)
%TRANSMIT_SPECTRUM Spectral density of the transmitted waveforms and checks of the estimator.
%   Streams: conventional OFDM and F-OFDM without a limiter, and each transmit chain at the
%   clipping ratio selected in selected_gamma.csv. The density is estimated from
%   p.psdFrames frames (seed p.calSeed) by averaged periodograms of Kaiser-windowed segments
%   (BAND_LEVELS, WELCH_SUM), each frame without its first and last p.Nb sampling periods.
%   Writes to resultsDir:
%   psd.csv        density of each stream for frequencies within half the sampling rate, in
%                  decibels relative to the mean density for |f| <= p.inBandHz
%   psd_check.csv  guard-band level of OFDM and F-OFDM without a limiter for the default
%                  estimator, a four times longer segment, a larger window shape parameter,
%                  and the closed-form density of the cyclic-prefixed block stream
%                  (ANALYTIC_GUARD_LEVEL)
p = sim_params();
ctx = sim_context(p);
cfg = configurations(p, resultsDir);
names = cfg.name;
chainOf = cfg.chain;
gamma = cfg.gamma;
nC = numel(names);
% Estimator variants of the check: segment length and shape parameter.
variants = [p.psdSegment, p.kaiserBeta; 4 * p.psdSegment, p.kaiserBeta; p.psdSegment, 40];
nV = size(variants, 1);
Tmax = p.B * p.Nb + p.L - 1;

acc = cell(1, nC);
count = zeros(1, nC);
accV = cell(2, nV);
countV = zeros(2, nV);
stream = RandStream('mt19937ar', 'Seed', p.calSeed);
for f = 1:p.psdFrames
    bits = randi(stream, [0 1], 4 * p.Nd, p.B);
    [s0, ~, ~, X] = ofdm_frame(p, bits);
    fr = {prepare_frame(p, ctx, s0, ones(p.B, 1), zeros(Tmax, 1), false), ...
        prepare_frame(p, ctx, s0, ones(p.B, 1), zeros(Tmax, 1), true)};
    for c = 1:nC
        ic = find(strcmp(p.chains, chainOf{c}));
        [~, xn, rate] = transmit_stream(p, ctx, fr{1}, fr{2}, chainOf{c}, gamma(c), X);
        T = numel(fr{p.chainFiltered(ic) + 1}.u);
        t = xn(rate * p.Nb + 1:rate * (T - p.Nb));
        if isempty(acc{c})
            acc{c} = zeros(rate * p.psdSegment, 1);
        end
        [acc{c}, count(c)] = welch_sum(acc{c}, count(c), t, kaiser_window(rate * p.psdSegment, p.kaiserBeta));
        if c <= 2
            for v = 1:nV
                if isempty(accV{c, v})
                    accV{c, v} = zeros(variants(v, 1), 1);
                end
                [accV{c, v}, countV(c, v)] = welch_sum(accV{c, v}, countV(c, v), t, ...
                    kaiser_window(variants(v, 1), variants(v, 2)));
            end
        end
    end
end

psd = zeros(p.psdSegment, nC);
for c = 1:nC
    rate = numel(acc{c}) / p.psdSegment;
    [psd(:, c), freq] = band_levels(p, acc{c} / count(c), rate);
end
fid = fopen(fullfile(resultsDir, 'psd.csv'), 'w');
fprintf(fid, 'freq_mhz');
for c = 1:nC
    fprintf(fid, ',%s_db', names{c});
end
fprintf(fid, '\n');
fprintf(fid, ['%.4f' repmat(',%.3f', 1, nC) '\n'], [freq / 1e6, round(1e4 * log10(psd)) / 1e3 + 0].');
fclose(fid);

fid = fopen(fullfile(resultsDir, 'psd_check.csv'), 'w');
fprintf(fid, 'stream,estimator,segment,kaiser_beta,frames,guard_band_db\n');
for c = 1:2
    for v = 1:nV
        q = p;
        q.psdSegment = variants(v, 1);
        [~, ~, guardDb] = band_levels(q, accV{c, v} / countV(c, v), 1);
        fprintf(fid, '%s,welch,%d,%d,%d,%s\n', names{c}, variants(v, 1), variants(v, 2), p.psdFrames, ...
            num_text(guardDb, 4));
    end
    if c == 1
        h = [];
    else
        h = ctx.h;
    end
    fprintf(fid, '%s,closed_form,NaN,NaN,NaN,%s\n', names{c}, num_text(analytic_guard_level(p, h), 4));
end
fclose(fid);
end

function guardDb = analytic_guard_level(p, h)
%ANALYTIC_GUARD_LEVEL Guard-band level of the closed-form density of the block stream.
%   For independent unit-power data symbols, the density of the cyclic-prefixed block stream
%   is proportional to the sum over the data subcarriers k of
%   |sin(pi*Nb*d)/sin(pi*d)|^2 with d = k/N - f/Fs, multiplied by |H(f)|^2 when the stream
%   passes the subband filter h. The level is the mean over p.guardHz <= |f| <= Fs/2
%   relative to the mean over |f| <= p.inBandHz on a grid of 2^16 frequencies.
nf = 2^16;
fr = (-nf / 2:nf / 2 - 1).' / nf;                  % frequency relative to the sampling rate
S = zeros(nf, 1);
for k = p.k.'
    d = k / p.N - fr;
    num = sin(pi * p.Nb * d);
    den = sin(pi * d);
    term = (num ./ den).^2;
    term(abs(den) < 1e-15) = p.Nb^2;
    S = S + term;
end
if ~isempty(h)
    H = fftshift(fft(h, nf));
    S = S .* abs(H).^2;
end
f = fr * p.Fs;
guardDb = 10 * log10(mean(S(abs(f) >= p.guardHz)) / mean(S(abs(f) <= p.inBandHz)));
end
