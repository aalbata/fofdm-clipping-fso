function characterize_transmitters(resultsDir)
%CHARACTERIZE_TRANSMITTERS Power, attenuation factor, peaks, and guard band of every transmit chain.
%   Uses p.calFrames frames of random bits from the seed p.calSeed, which no run of the
%   bit-error simulation uses. For every chain and clipping ratio:
%   - average power of the transmitted waveform over the counted blocks;
%   - attenuation factor: real part of the correlation of the noise-free FFT outputs
%     (unit channel gain, divided by the response of the filter cascade and the timing
%     advance) with the data symbols, divided by the symbol energy, and the in-band
%     signal-to-distortion ratio of the outputs relative to that factor;
%   - peak power: the transmitted waveform is reconstructed at p.peakRate times the sampling
%     rate by band-limited interpolation, and the instantaneous power exceeded by a fraction
%     q of its samples over the counted blocks is reported relative to the average power
%     for each q in p.peakProbs, with a 95% bootstrap interval over frames, together with
%     the largest sample;
%   - guard-band level and power beyond half the sampling rate (BAND_LEVELS), from the
%     frames without their first and last p.Nb sampling periods.
%   Writes to resultsDir:
%   calibration.csv       one row per chain and clipping ratio
%   peak_convergence.csv  the peak powers at p.peakRate and p.peakRateCheck times the
%                         sampling rate from the first p.checkFrames frames, for the clipping
%                         ratios p.checkGamma
%   ccdf.csv              fraction of samples whose power exceeds the average power by more
%                         than a given level, for every chain and clipping ratio
p = sim_params();
ctx = sim_context(p);
nC = numel(p.chains);
nG = numel(p.gammaGrid);
nF = p.calFrames;
nQ = numel(p.peakProbs);
if nQ ~= 3 || p.peakProbs(2) ~= p.peakProb
    error('characterize_transmitters:peakProb', 'p.peakProbs must hold three values with p.peakProb second.');
end
edges = (-30:p.levelBinDb:18).';                   % level relative to the limiter input power (dB)
nE = numel(edges) - 1;
Tmax = p.B * p.Nb + p.L - 1;

alphaNum = zeros(nC, nG);
zEnergy = zeros(nC, nG);
xEnergy = zeros(nC, nG);
hist16 = zeros(nC, nG, nF, nE);
pSum = zeros(nC, nG, nF);
pCount = zeros(nC, nG, nF);
pMax = zeros(nC, nG);
histCheck = zeros(nC, nG, 2, nE);                  % rate index 1: p.peakRate, 2: p.peakRateCheck
pSumCheck = zeros(nC, nG, 2);
pCountCheck = zeros(nC, nG, 2);
acc = cell(nC, nG);
count = zeros(nC, nG);
win = {kaiser_window(p.psdSegment, p.kaiserBeta), kaiser_window(p.os * p.psdSegment, p.kaiserBeta)};
rates = [p.peakRate, p.peakRateCheck];
isCheck = ismember(round(10 * p.gammaGrid), round(10 * p.checkGamma));

stream = RandStream('mt19937ar', 'Seed', p.calSeed);
tStart = tic;
for f = 1:nF
    bits = randi(stream, [0 1], 4 * p.Nd, p.B);
    [s0, ~, ~, X] = ofdm_frame(p, bits);
    Xc = X(:, 2:p.M + 1);
    fr = {prepare_frame(p, ctx, s0, ones(p.B, 1), zeros(Tmax, 1), false), ...
        prepare_frame(p, ctx, s0, ones(p.B, 1), zeros(Tmax, 1), true)};
    for ic = 1:nC
        frc = fr{p.chainFiltered(ic) + 1};
        s2 = ctx.chainSigma2(ic);
        T = numel(frc.u);
        for ig = 1:nG
            [xs, xn, rate] = transmit_stream(p, ctx, fr{1}, fr{2}, p.chains{ic}, p.gammaGrid(ig), X);
            Z = equalized_symbols(p, ctx, frc, xs);
            alphaNum(ic, ig) = alphaNum(ic, ig) + real(sum(Z(:) .* conj(Xc(:))));
            zEnergy(ic, ig) = zEnergy(ic, ig) + sum(abs(Z(:)).^2);
            xEnergy(ic, ig) = xEnergy(ic, ig) + sum(abs(Xc(:)).^2);

            % Peaks of the waveform reconstructed at p.peakRate times the sampling rate.
            [h, ps, pc, pm] = power_histogram(p, xn, p.peakRate, frc.D, s2, edges);
            hist16(ic, ig, f, :) = h;
            pSum(ic, ig, f) = ps;
            pCount(ic, ig, f) = pc;
            pMax(ic, ig) = max(pMax(ic, ig), pm);
            if f <= p.checkFrames && isCheck(ig)
                histCheck(ic, ig, 1, :) = reshape(histCheck(ic, ig, 1, :), 1, []) + h(:).';
                pSumCheck(ic, ig, 1) = pSumCheck(ic, ig, 1) + ps;
                pCountCheck(ic, ig, 1) = pCountCheck(ic, ig, 1) + pc;
                [hc, psc, pcc] = power_histogram(p, xn, rates(2), frc.D, s2, edges);
                histCheck(ic, ig, 2, :) = reshape(histCheck(ic, ig, 2, :), 1, []) + hc(:).';
                pSumCheck(ic, ig, 2) = pSumCheck(ic, ig, 2) + psc;
                pCountCheck(ic, ig, 2) = pCountCheck(ic, ig, 2) + pcc;
            end

            % Spectral density of the frame without its first and last p.Nb sampling periods.
            if isempty(acc{ic, ig})
                acc{ic, ig} = zeros(rate * p.psdSegment, 1);
            end
            t = xn(rate * p.Nb + 1:rate * (T - p.Nb));
            [acc{ic, ig}, count(ic, ig)] = welch_sum(acc{ic, ig}, count(ic, ig), t, win{1 + (rate > 1)});
        end
    end
    fprintf('characterization frame %d/%d (%.0f s)\n', f, nF, toc(tStart));
end

% Bootstrap over frames for the peak powers.
boot = RandStream('mt19937ar', 'Seed', p.bootSeed);
draws = randi(boot, nF, nF, p.bootDraws);
qName = arrayfun(@(q) sprintf('papr_1e%d', round(-log10(q))), p.peakProbs, 'UniformOutput', false);
fid = fopen(fullfile(resultsDir, 'calibration.csv'), 'w');
fprintf(fid, 'chain,gamma,limiter_sigma2,tx_power,alpha,alpha_gauss,power_gauss,inband_sdr_db');
fprintf(fid, ',%s_db', qName{:});
fprintf(fid, ',papr_max_db');
for iq = 1:nQ
    fprintf(fid, ',%s_lo_db,%s_hi_db', qName{iq}, qName{iq});
end
fprintf(fid, ',guard_band_db,outside_band_db\n');
cfid = fopen(fullfile(resultsDir, 'ccdf.csv'), 'w');
fprintf(cfid, 'chain,gamma,level_db,exceedance\n');
for ic = 1:nC
    s2 = ctx.chainSigma2(ic);
    for ig = 1:nG
        gamma = p.gammaGrid(ig);
        alpha = alphaNum(ic, ig) / xEnergy(ic, ig);
        sdr = alpha^2 * xEnergy(ic, ig) / (zEnergy(ic, ig) - alpha^2 * xEnergy(ic, ig));
        H = reshape(sum(hist16(ic, ig, :, :), 3), [], 1);
        total = sum(pCount(ic, ig, :));
        txPower = sum(pSum(ic, ig, :)) / total;
        offset = 10 * log10(txPower / s2);         % average power relative to the limiter input (dB)
        papr = exceedance_level(H, edges, total, p.peakProbs) - offset;
        bootPapr = zeros(p.bootDraws, nQ);
        for ib = 1:p.bootDraws
            fs = draws(:, ib);
            Hb = reshape(sum(hist16(ic, ig, fs, :), 3), [], 1);
            tb = sum(pCount(ic, ig, fs));
            bootPapr(ib, :) = exceedance_level(Hb, edges, tb, p.peakProbs) - 10 * log10(sum(pSum(ic, ig, fs)) / tb / s2);
        end
        bootPapr = sort(bootPapr, 1);
        lo = bootPapr(round(0.025 * p.bootDraws), :);
        hi = bootPapr(round(0.975 * p.bootDraws), :);
        rate = 1 + (numel(acc{ic, ig}) > p.psdSegment) * (p.os - 1);
        [~, ~, guardDb, outsideDb] = band_levels(p, acc{ic, ig} / count(ic, ig), rate);
        ref = clipping_reference(gamma);
        if strcmp(p.chains{ic}, 'rcf') && ~isinf(gamma)
            alphaGauss = NaN;
            powerGauss = NaN;
        elseif strcmp(p.chains{ic}, 'clip_filter') && ~isinf(gamma)
            alphaGauss = ref.alpha;                % the filter after the limiter is linear
            powerGauss = NaN;                      % the filter removes part of the distortion
        elseif strcmp(p.chains{ic}, 'clip_filter')
            alphaGauss = ref.alpha;
            powerGauss = ctx.sigma2(2);            % without a limiter the chain transmits F-OFDM
        else
            alphaGauss = ref.alpha;
            powerGauss = s2 * ref.power;
        end
        fprintf(fid, '%s,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%s', p.chains{ic}, gamma_text(gamma), s2, txPower, ...
            alpha, alphaGauss, powerGauss, num_text(10 * log10(sdr), 4));
        fprintf(fid, ',%s', num_text(papr(1), 4), num_text(papr(2), 4), num_text(papr(3), 4), ...
            num_text(10 * log10(pMax(ic, ig) / txPower), 4));
        for iq = 1:nQ
            fprintf(fid, ',%s,%s', num_text(lo(iq), 4), num_text(hi(iq), 4));
        end
        fprintf(fid, ',%s,%s\n', num_text(guardDb, 4), num_text(outsideDb, 4));

        % Exceedance of levels relative to the average power, on a grid of 0.05 dB, interpolated
        % linearly between the histogram edges.
        above = [flipud(cumsum(flipud(H))); 0];
        for lev = 0:0.05:14
            x = (lev + offset - edges(1)) / p.levelBinDb;
            j = floor(x) + 1;
            if j < 1 || j > nE
                continue
            end
            n = above(j) - (x - (j - 1)) * (above(j) - above(j + 1));
            if n > 0
                fprintf(cfid, '%s,%s,%.2f,%.4e\n', p.chains{ic}, gamma_text(gamma), lev, n / total);
            end
        end
    end
end
fclose(fid);
fclose(cfid);

fid = fopen(fullfile(resultsDir, 'peak_convergence.csv'), 'w');
fprintf(fid, 'chain,gamma,rate');
fprintf(fid, ',%s_db', qName{:});
fprintf(fid, '\n');
for ic = 1:nC
    for ig = find(isCheck)
        for ir = 1:2
            total = pCountCheck(ic, ig, ir);
            offset = 10 * log10(pSumCheck(ic, ig, ir) / total / ctx.chainSigma2(ic));
            papr = exceedance_level(reshape(histCheck(ic, ig, ir, :), [], 1), edges, total, p.peakProbs) - offset;
            fprintf(fid, '%s,%s,%d,%s,%s,%s\n', p.chains{ic}, gamma_text(p.gammaGrid(ig)), rates(ir), ...
                num_text(papr(1), 4), num_text(papr(2), 4), num_text(papr(3), 4));
        end
    end
end
fclose(fid);
end
