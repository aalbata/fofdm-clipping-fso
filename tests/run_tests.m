function run_tests()
%RUN_TESTS Checks of the building blocks of the simulation and of the analysis.
%   Prints one line per check and stops with an error at the first failure.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
p = sim_params();
ctx = sim_context(p);
check(p.os == 4 && p.nPad == 2^19 && p.M == 60 && p.delta == 205 && p.rcfIterations == 4 && ...
    p.peakRate == 16 && p.peakProb == 1e-4, 'simulation constants');

% Chebyshev window: symmetry, peak, and equal sidelobes at the specified level.
w = dolph_chebyshev(p.L, p.sidelobeDb);
check(abs(max(w) - 1) < 1e-12 && max(abs(w - flipud(w))) < 1e-9, 'window is symmetric with unit maximum');
W = 20 * log10(abs(fft(w, 64 * p.L)) / sum(w));
firstNull = find(diff(W(1:32 * p.L)) > 0, 1);
check(abs(max(W(firstNull:32 * p.L)) + p.sidelobeDb) < 0.05, 'window sidelobes are 75 dB below the main lobe');
if exist('chebwin', 'file') == 2
    check(max(abs(w - chebwin(p.L, p.sidelobeDb))) < 1e-9, 'window equals chebwin of the Signal Processing Toolbox');
end

% Subband filter: symmetry, gain, passband, half-amplitude point, and stopband.
Hf = abs(fft(ctx.h, p.N));
check(abs(sum(ctx.h) - 1) < 1e-12 && max(abs(ctx.h - flipud(ctx.h))) < 1e-12, ...
    'filter is symmetric with unit gain at zero frequency');
check(max(abs(20 * log10(Hf(1:1501)))) < 0.01, 'filter passband is flat within 0.01 dB for |k| <= 1500');
check(abs(Hf(1540) - 0.5) < 1e-3, 'filter has half amplitude at the band edge plus the tone offset');
check(abs(20 * log10(Hf(1537)) + 3.657) < 0.01, 'filter attenuates the outermost data subcarrier by 3.66 dB');
check(max(20 * log10(Hf(1601:p.N / 2 + 1))) < -80, 'filter response is below -80 dB for |k| >= 1600');
check(abs(ctx.sigma2(1) - 0.75) < 1e-15 && abs(ctx.sigma2(2) - 0.7482979198) < 1e-9, ...
    'stream power without and with the filter');
check(isequal(ctx.chainSigma2, ctx.sigma2([1 2 1 1])), 'limiter input power of each chain');

% Gaussian clipping reference against numerical integration over the Rayleigh envelope and
% against complex Gaussian samples.
for gamma = [0.7 1.5 3]
    ref = clipping_reference(gamma);
    alphaInt = integral(@(r) min(r, gamma) .* r .* 2 .* r .* exp(-r.^2), 0, Inf, 'Waypoints', gamma, ...
        'AbsTol', 1e-14, 'RelTol', 1e-12);
    powerInt = integral(@(r) min(r, gamma).^2 .* 2 .* r .* exp(-r.^2), 0, Inf, 'Waypoints', gamma, ...
        'AbsTol', 1e-14, 'RelTol', 1e-12);
    check(abs(alphaInt - ref.alpha) < 1e-10 && abs(powerInt - ref.power) < 1e-10, ...
        sprintf('closed-form attenuation factor and output power at gamma = %.1f', gamma));
    check(abs(10 * log10(ref.papr) - 10 * log10(gamma^2 / powerInt)) < 1e-9, sprintf('PAPR at gamma = %.1f', gamma));
end
stream = RandStream('mt19937ar', 'Seed', 1);
u = (randn(stream, 2e6, 1) + 1j * randn(stream, 2e6, 1)) / sqrt(2);
for gamma = [0.7 1.5 3]
    ref = clipping_reference(gamma);
    c = envelope_limiter(u, gamma);
    d = c - ref.alpha * u;
    check(abs(real(mean(c .* conj(u))) - ref.alpha) < 3e-3, sprintf('attenuation factor at gamma = %.1f', gamma));
    check(abs(mean(abs(c).^2) - ref.power) < 3e-3, sprintf('output power at gamma = %.1f', gamma));
    check(abs(mean(abs(d).^2) - ref.distortion) < 3e-3, sprintf('distortion power at gamma = %.1f', gamma));
    check(max(abs(c)) <= gamma * (1 + 1e-12), sprintf('output magnitude bounded at gamma = %.1f', gamma));
end
ref = clipping_reference(Inf);
check(ref.alpha == 1 && ref.power == 1 && ref.distortion == 0, 'reference without limiter');

% Exceedance level from a histogram: unit-mean exponential power, P(power > x) = exp(-x).
pw = abs(u).^2;
edges = (0:0.01:18).';
lev = exceedance_level(histcounts(10 * log10(pw), edges), edges, numel(pw), [1e-2 1e-3]);
check(all(abs(lev - 10 * log10(-log([1e-2 1e-3]))) < [0.01 0.03]), 'exceedance level of an exponential power');

% Rate change: a tone inside the band is interpolated exactly, and decimation inverts interpolation.
n = 4096;
tone = exp(2j * pi * 37 * (0:n - 1).' / n);
check(max(abs(resample_stream(tone, 4 * n) - exp(2j * pi * 37 * (0:4 * n - 1).' / (4 * n)))) < 1e-10, ...
    'interpolation of a tone');
x = u(1:n);
check(max(abs(resample_stream(resample_stream(x, 4 * n), n) - x)) < 1e-10, 'decimation inverts interpolation');
Xb = fft(x);
Xb(n / 4 + 1:3 * n / 4) = 0;                       % band-limited stream: upper half of the band removed
xb = ifft(Xb);
xo = resample_stream(xb, 4 * n);
check(abs(mean(abs(xo).^2) - mean(abs(xb).^2)) < 1e-10, 'interpolation preserves the power of a band-limited stream');
check(max(abs(xo(1:4:end) - xb)) < 1e-10, 'interpolation keeps the original samples');
high = exp(2j * pi * (3 * n) * (0:4 * n - 1).' / (4 * n));         % tone beyond half the lower sampling rate
check(max(abs(resample_stream(high, n))) < 1e-10, 'decimation removes a tone outside the lower band');

% Bit-error counting from crossing scales against direct decisions at each noise scale.
nSym = 20000;
lv = randi(stream, [0 3], nSym, 2);
Xt = ((2 * lv(:, 1) - 3) + 1j * (2 * lv(:, 2) - 3)) / sqrt(10);
Zs = 0.9 * Xt + 0.05 * (randn(stream, nSym, 1) + 1j * randn(stream, nSym, 1));
Zn = (randn(stream, nSym, 1) + 1j * randn(stream, nSym, 1)) / sqrt(2);
Zn(1:10) = 0;
scale = [0, 10.^(-(0:0.25:30) / 20)];
check(isequal(count_bit_errors(Zs, Zn, scale, lv(:, 1), lv(:, 2)), direct_bit_errors(Zs, Zn, scale, lv(:, 1), lv(:, 2))), ...
    'bit-error counting equals direct decisions at every noise scale');
% Exact ties: at unit noise scale the noise cancels every component.
Zt = Zs(1:2000);
lt = lv(1:2000, :);
st = [0.5 1 2];
check(isequal(count_bit_errors(Zt, -Zt, st, lt(:, 1), lt(:, 2)), direct_bit_errors(Zt, -Zt, st, lt(:, 1), lt(:, 2))), ...
    'bit-error counting at exact sign ties');

% Closed-form 16-QAM bit error rate and required SNR.
check(abs(qam16_ber(0) - 0.5) < 1e-15, 'closed-form BER at zero SNR');
snr = 0:5;
ber = 10.^(-(0:5) / 2);
check(abs(required_snr(snr, ber, 10^(-1.25)) - 2.5) < 1e-9, 'required SNR on a straight line of log BER');
check(isnan(required_snr(snr, ber, 1e-9)), 'required SNR is NaN when the target is not reached');
check(isnan(required_snr(snr, [1 1 1 0 0 0], 1e-2)), 'required SNR is NaN when the BER drops to zero at the crossing');
check(abs(required_snr(0:2, [0 1e-2 1e-4], 1e-3) - 1.5) < 1e-9, ...
    'required SNR falls back to the straight line without a third usable point');
curve = qam16_ber(10.^(p.snrGrid / 10) * p.N / p.Nd);
for target = p.berTargets
    exact = fzero(@(s) log10(qam16_ber(10^(s / 10) * p.N / p.Nd)) - log10(target), [0 30]);
    check(abs(required_snr(p.snrGrid, curve, target) - exact) < 5e-4, ...
        sprintf('required SNR of the closed form at BER %g within 0.0005 dB', target));
    check(abs(required_snr(p.snrGrid, curve, target, 'linear') - exact) < 0.01, ...
        sprintf('linear interpolation of the closed form at BER %g within 0.01 dB', target));
end

% Required-SNR table: curves shifted by a known amount per clipping ratio.
nG = numel(p.gammaGrid);
nC = numel(p.chains);
shift = 0.1 * (1:nG);
S.ber = zeros(numel(p.records), nC, nG, numel(p.snrGrid));
for ig = 1:nG
    S.ber(:, :, ig, :) = repmat(reshape(qam16_ber(10.^((p.snrGrid - shift(ig)) / 10) * p.N / p.Nd), 1, 1, 1, []), ...
        numel(p.records), nC);
end
S.berSeed = repmat(S.ber, 1, 1, 1, 1, numel(p.seeds));
R = required_table(p, S);
exact = fzero(@(s) log10(qam16_ber(10^(s / 10) * p.N / p.Nd)) - log10(p.berTargets(1)), [0 30]);
check(max(abs(reshape(R.req(1, 2, :, 1), 1, []) - exact - shift)) < 5e-4, 'required-SNR table follows the shift of the curves');
check(max(abs(R.reqSeed(:) - repmat(R.req(:), numel(p.seeds), 1))) < 1e-12, 'required-SNR table per run');
check(max(abs(R.req1dB(:) - R.req(:))) < 0.005, 'required SNR from the whole-decibel grid');

% One frame through every transmit chain and the receiver.
Tmax = p.B * p.Nb + p.L - 1;
bits = randi(stream, [0 1], 4 * p.Nd, p.B);
w0 = (randn(stream, Tmax, 1) + 1j * randn(stream, Tmax, 1)) / sqrt(2);
ab = 0.8 + 0.4 * rand(stream, p.B, 1);
[s0, iq, qq, X] = ofdm_frame(p, bits);
iq = iq(:, 2:p.M + 1);
qq = qq(:, 2:p.M + 1);
Xc = X(:, 2:p.M + 1);
check(abs(mean(abs(Xc(:)).^2) - 1) < 5e-3, 'transmitted symbols have unit average power');
nBits = 4 * p.Nd * p.M;
fr = {prepare_frame(p, ctx, s0, ab, w0, false), prepare_frame(p, ctx, s0, ab, w0, true)};
check(numel(fr{1}.uo) == p.os * p.nPad && numel(fr{2}.uo) == p.os * p.nPad, 'interpolated frame length');
minSdr = [200 70 70 200];                          % without filter: rounding only; with filter: block interference
for ic = 1:nC
    chain = p.chains{ic};
    frc = fr{p.chainFiltered(ic) + 1};
    [xs, xn, rate] = transmit_stream(p, ctx, fr{1}, fr{2}, chain, Inf, X);
    check(rate == 1 && isequal(xs, frc.u) && numel(xn) == p.nPad, ['stream without a limiter, ' chain]);
    Z = equalized_symbols(p, ctx, frc, xs);
    check(count_bit_errors(Z, frc.Zn, 0, iq, qq) == 0, ['no bit errors without noise and limiter, ' chain]);
    sdr = 10 * log10(sum(abs(Xc(:)).^2) / sum(abs(Z(:) - Xc(:)).^2));
    check(sdr > minSdr(ic), ['equalized symbols equal the transmitted symbols, ' chain]);
    for gamma = [1.5 3]
        ref = clipping_reference(gamma);
        A = gamma * sqrt(ctx.chainSigma2(ic));
        [xs, xn, rate] = transmit_stream(p, ctx, fr{1}, fr{2}, chain, gamma, X);
        Z = equalized_symbols(p, ctx, frc, xs);
        alphaMeas = real(sum(Z(:) .* conj(Xc(:)))) / sum(abs(Xc(:)).^2);
        tag = sprintf('gamma = %.1f, %s', gamma, chain);
        switch chain
            case {'clip', 'filter_clip'}
                check(rate == p.os && max(abs(xn)) <= A * (1 + 1e-12), ['limiter output bounded, ' tag]);
                check(abs(alphaMeas - ref.alpha) < 3e-3, ['attenuation factor equals the Gaussian reference, ' tag]);
                pOut = mean(abs(xn(p.os * (frc.D + p.Nb) + (1:p.os * p.M * p.Nb))).^2);
                check(abs(pOut / (ctx.chainSigma2(ic) * ref.power) - 1) < 5e-3, ['limiter output power, ' tag]);
            case 'clip_filter'
                check(rate == 1 && numel(xs) == numel(frc.u), ['filtered stream length, ' tag]);
                check(abs(alphaMeas - ref.alpha) < 3e-3, ['attenuation factor equals the Gaussian reference, ' tag]);
            case 'rcf'
                blocks = reshape(xs, p.Nb, p.B);
                F = fft(blocks(p.Ncp + 1:end, :), [], 1);
                outside = true(p.N, 1);
                outside(p.bins) = false;
                check(max(abs(F(outside, :)), [], 'all') < 1e-9 * max(abs(F(:))), ...
                    ['no components outside the data subcarriers, ' tag]);
                check(isequal(blocks(1:p.Ncp, :), blocks(p.N + 1:end, :)), ['cyclic prefix, ' tag]);
                xo = resample_stream(s0, p.os * numel(s0));
                xr = resample_stream(xs, p.os * numel(xs));
                check(max(abs(xr)) < max(abs(xo)), ['repeated clipping and filtering lowers the peak, ' tag]);
        end

        % Peak-limited model simulated directly: waveform scaled to the peak power limit, noise
        % added, complete receiver. The result equals the average-power result at an SNR
        % lower by the PAPR.
        pAvg = mean(abs(xs(frc.D + p.Nb + (1:p.M * p.Nb))).^2);
        pPeak = 3 * pAvg;                          % any peak power gives the same identity
        psnr = 16;
        snrEq = psnr - 10 * log10(pPeak / pAvg);
        y = frc.an .* xs / sqrt(pPeak) + 10^(-psnr / 20) * w0(1:numel(xs));
        if frc.filt
            y = filter_stream(y, ctx.g);
        end
        Zhat = receive_frame(p, y, frc.Dc) ./ (frc.eqz * alphaMeas / sqrt(pPeak));
        eDirect = count_bit_errors(Zhat, zeros(size(Zhat)), 0, iq, qq);
        eShift = count_bit_errors(Z / alphaMeas, frc.Zn, 10^(-snrEq / 20) * sqrt(pAvg) / alphaMeas, iq, qq);
        check(eDirect > 100 && abs(eDirect - eShift) <= 2, ['peak-limited model equals the shifted average-power model, ' tag]);
    end
end

% Unfiltered OFDM without limiter and with unit gains against the closed form.
fr1 = prepare_frame(p, ctx, s0, ones(p.B, 1), w0, false);
Z = equalized_symbols(p, ctx, fr1, fr1.u);
snr = [8 12 16];
e = count_bit_errors(Z, fr1.Zn, 10.^(-snr / 20) * sqrt(p.Nd / p.N), iq, qq);
expected = nBits * qam16_ber(10.^(snr / 10) * p.N / p.Nd);
check(all(abs(e - expected) < 4 * sqrt(expected)), 'BER in white Gaussian noise matches the closed form');
fprintf('All checks passed.\n');
end

function e = direct_bit_errors(Zs, Zn, scale, iq, qq)
%DIRECT_BIT_ERRORS Bit errors of hard 16-QAM decisions evaluated separately at each noise scale.
s = sqrt(10) / 2;
v0 = [real(Zs(:)); imag(Zs(:))] * s;
n0 = [real(Zn(:)); imag(Zn(:))] * s;
level = [iq(:); qq(:)];
sentPositive = level >= 2;
sentInner = level == 1 | level == 2;
e = zeros(1, numel(scale));
for i = 1:numel(scale)
    v = v0 + scale(i) * n0;
    e(i) = nnz(xor(v > 0, sentPositive)) + nnz(xor(abs(v) < 1, sentInner));
end
end

function check(condition, name)
%CHECK Report a passed check or stop at a failed one.
if ~condition
    error('run_tests:failed', 'FAILED: %s', name);
end
fprintf('passed: %s\n', name);
end
