function implementation_checks(dataDir, resultsDir)
%IMPLEMENTATION_CHECKS Checks of simplifications of the transmitter and receiver models.
%   Writes implementation_checks.csv to resultsDir with one row per check and configuration:
%   gain_order    the simulation multiplies the waveform within the sampling bandwidth by the
%                 block gains at the sampling rate. The check multiplies the transmitted
%                 waveform at p.os times the sampling rate by the block gains and restricts
%                 the product to the sampling bandwidth afterwards. The value is the ratio in
%                 decibels of the energy of the noise-free FFT outputs to the energy of their
%                 difference between the two orders, for p.B gains of the rainy record taken
%                 at every numel(p.seeds)-th sample from the first (the subsequence used by
%                 run 1), at the selected clipping ratio of each chain.
%   noise_power   mean power of the receiver noise at the FFT output on the data
%                 subcarriers, multiplied by the squared magnitude of the known response and
%                 divided by the response of the filter cascade, for unit-variance noise and
%                 unit gain; the expected value is one.
%   hard_limit    the peak-limited model passes the fraction p.peakProb of the samples above
%                 the peak power unchanged. The check limits the envelope of the transmitted
%                 waveform reconstructed at p.peakRate times the sampling rate to the square
%                 root of the peak power in calibration.csv, restricts it to the sampling
%                 bandwidth, and reports the ratio in decibels of the energy of the noise-free
%                 FFT outputs to the energy of the change caused by the limiting, over the first
%                 p.checkFrames characterization frames, for the six configurations.
%   clipped_fraction  fraction of the reconstructed samples limited in the hard_limit check.
%   lowpass_papr  PAPR (peak power at p.peakProb, reconstruction at p.peakRate times the
%                 sampling rate) of the clip chain after its output is restricted to the
%                 sampling bandwidth, without a subband filter, over the p.calFrames
%                 characterization frames, at the selected clipping ratios of clip and
%                 clip_filter. Together with calibration.csv it separates the peak regrowth
%                 caused by the bandwidth restriction from that caused by the subband filter.
p = sim_params();
ctx = sim_context(p);
cfg = configurations(p, resultsDir);
[~, ~, cal] = read_calibration(p, resultsDir);
Tmax = p.B * p.Nb + p.L - 1;
fid = fopen(fullfile(resultsDir, 'implementation_checks.csv'), 'w');
fprintf(fid, 'check,configuration,gamma,value\n');

% Order of the block gains and the restriction to the sampling bandwidth; noise power.
rec = load_power_record(dataDir, 'rain');
ab = rec.a(1:numel(p.seeds):numel(p.seeds) * p.B);
stream = RandStream('mt19937ar', 'Seed', 77);
bits = randi(stream, [0 1], 4 * p.Nd, p.B);
w0 = (randn(stream, Tmax, 1) + 1j * randn(stream, Tmax, 1)) / sqrt(2);
[s0, ~, ~, X] = ofdm_frame(p, bits);
fr = {prepare_frame(p, ctx, s0, ab, w0, false), prepare_frame(p, ctx, s0, ab, w0, true)};
fr1 = {prepare_frame(p, ctx, s0, ones(p.B, 1), w0, false), prepare_frame(p, ctx, s0, ones(p.B, 1), w0, true)};
for ic = 1:numel(p.chains)
    chain = p.chains{ic};
    c = find(strcmp(cfg.name, chain));
    gamma = cfg.gamma(c);
    frc = fr{p.chainFiltered(ic) + 1};
    T = numel(frc.u);
    [xs, xn, rate] = transmit_stream(p, ctx, fr{1}, fr{2}, chain, gamma, X);
    Z = equalized_symbols(p, ctx, frc, xs);
    if rate == 1
        xn = resample_stream(xn, p.os * p.nPad);
    end
    an4 = repelem([frc.an; frc.an(end) * ones(p.nPad - T, 1)], p.os);
    y = resample_stream(an4 .* xn, p.nPad);
    y = y(1:T);
    if frc.filt
        y = filter_stream(y, ctx.g);
    end
    Zalt = receive_frame(p, y, frc.Dc) ./ frc.eqz;
    fprintf(fid, 'gain_order,%s,%s,%.2f\n', chain, gamma_text(gamma), ...
        10 * log10(sum(abs(Z(:)).^2) / sum(abs(Z(:) - Zalt(:)).^2)));
    f1 = fr1{p.chainFiltered(ic) + 1};
    Gk = ones(p.Nd, 1);
    if f1.filt
        Gk = ctx.G;
    end
    fprintf(fid, 'noise_power,%s,Inf,%.4f\n', chain, mean(abs(f1.Zn).^2 .* abs(f1.eqz).^2 ./ Gk, 'all'));
end

% Hard limit at the peak power, and the clip chain restricted to the sampling bandwidth.
nCfg = numel(cfg.name);
zE = zeros(1, nCfg);
dE = zeros(1, nCfg);
nClip = zeros(1, nCfg);
nAll = zeros(1, nCfg);
lowGamma = cfg.gamma(strcmp(cfg.name, 'clip') | strcmp(cfg.name, 'clip_filter'));
lowGamma = unique(lowGamma);
edges = (-30:p.levelBinDb:18).';
H = zeros(numel(edges) - 1, numel(lowGamma));
pS = zeros(1, numel(lowGamma));
pC = zeros(1, numel(lowGamma));
s2 = ctx.chainSigma2(strcmp(p.chains, 'clip'));
stream = RandStream('mt19937ar', 'Seed', p.calSeed);
for f = 1:p.calFrames
    bits = randi(stream, [0 1], 4 * p.Nd, p.B);
    [s0, ~, ~, X] = ofdm_frame(p, bits);
    fr = {prepare_frame(p, ctx, s0, ones(p.B, 1), zeros(Tmax, 1), false), ...
        prepare_frame(p, ctx, s0, ones(p.B, 1), zeros(Tmax, 1), true)};
    for k = 1:numel(lowGamma)
        xs = transmit_stream(p, ctx, fr{1}, fr{2}, 'clip', lowGamma(k), X);
        [h, ps, pc] = power_histogram(p, [xs; zeros(p.nPad - numel(xs), 1)], p.peakRate, 0, s2, edges);
        H(:, k) = H(:, k) + h(:);
        pS(k) = pS(k) + ps;
        pC(k) = pC(k) + pc;
    end
    if f > p.checkFrames
        continue
    end
    for c = 1:nCfg
        ic = cfg.ic(c);
        frc = fr{p.chainFiltered(ic) + 1};
        T = numel(frc.u);
        [xs, xn] = transmit_stream(p, ctx, fr{1}, fr{2}, cfg.chain{c}, cfg.gamma(c), X);
        row = cal.chain == cfg.chain{c} & (round(10 * cal.gamma) == round(10 * cfg.gamma(c)));
        aPk = sqrt(cal.tx_power(row) * 10^(cal.papr_1e4_db(row) / 10));
        xp = resample_stream(xn, p.peakRate * p.nPad);
        seg = p.peakRate * (frc.D + p.Nb) + (1:p.peakRate * p.M * p.Nb);
        nClip(c) = nClip(c) + nnz(abs(xp(seg)) > aPk);
        nAll(c) = nAll(c) + numel(seg);
        xl = resample_stream(envelope_limiter(xp, aPk), p.nPad);
        Z = equalized_symbols(p, ctx, frc, xs);
        Zl = equalized_symbols(p, ctx, frc, xl(1:T));
        zE(c) = zE(c) + sum(abs(Z(:)).^2);
        dE(c) = dE(c) + sum(abs(Zl(:) - Z(:)).^2);
    end
end
for c = 1:nCfg
    fprintf(fid, 'hard_limit,%s,%s,%.2f\n', cfg.name{c}, gamma_text(cfg.gamma(c)), 10 * log10(zE(c) / dE(c)));
    fprintf(fid, 'clipped_fraction,%s,%s,%.3e\n', cfg.name{c}, gamma_text(cfg.gamma(c)), nClip(c) / nAll(c));
end
for k = 1:numel(lowGamma)
    papr = exceedance_level(H(:, k), edges, pC(k), p.peakProb) - 10 * log10(pS(k) / pC(k) / s2);
    fprintf(fid, 'lowpass_papr,clip,%s,%.4f\n', gamma_text(lowGamma(k)), papr);
end
fclose(fid);
end
