function run_simulation(dataDir, resultsDir)
%RUN_SIMULATION Bit-error simulation of all configurations over both power records.
%   Reads the power records from dataDir and writes ber.csv, waveform_stats.csv, and
%   channel_stats.csv to resultsDir. The transmit power is normalized to one, so the noise
%   variance for an SNR of x dB is 10^(-x/10). Every configuration of a frame uses the same
%   bits, channel gains, and noise realization. Run i of a record uses the record samples
%   i, i + nSeed, i + 2*nSeed, ..., one sample per counted block, so that the runs together
%   use every sample once.
p = sim_params();
ctx = sim_context(p);
nR = numel(p.records);
nG = numel(p.gammaGrid);
nS = numel(p.snrGrid);
nSeed = numel(p.seeds);
noiseStd = 10.^(-p.snrGrid / 20);                  % noise standard deviation for unit transmit power
Tmax = p.B * p.Nb + p.L - 1;

errors = zeros(nR, 2, nG, nS, nSeed);              % record, filter, gamma, SNR, run
bitCount = zeros(nR, nSeed);
alphaNum = zeros(nR, nSeed, 2, nG);                % numerator of the measured attenuation factor
pSum = zeros(nR, nSeed, 2, nG);                    % sum over frames of the limiter output power
eX = zeros(nR, nSeed, 2, nG);                      % energy of the transmitted symbols
eErr = zeros(nR, nSeed, 2, nG);                    % energy of the noise-free symbol error
nFramesOf = zeros(nR, nSeed);
tStart = tic;

records = cell(1, nR);
for r = 1:nR
    records{r} = load_power_record(dataDir, p.records{r});
    if mod(numel(records{r}.a), nSeed * p.M) ~= 0
        error('run_simulation:recordLength', 'The length of record %s is not a multiple of %d.', ...
            p.records{r}, nSeed * p.M);
    end
end

chanFile = open_output(fullfile(resultsDir, 'channel_stats.csv'));
fprintf(chanFile, ['record,n_samples,sampling_time_s,duration_min,power_min_dbm,power_max_dbm,' ...
    'power_mean_dbm,gain_min_db,gain_max_db,scintillation_index\n']);

for r = 1:nR
    rec = records{r};
    Ts = rec.time(2) - rec.time(1);
    fprintf(chanFile, '%s,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', p.records{r}, ...
        numel(rec.power), Ts, (numel(rec.power) - 1) * Ts / 60, min(rec.power), max(rec.power), ...
        10 * log10(mean(10.^(rec.power / 10))), 10 * log10(min(rec.g)), 10 * log10(max(rec.g)), ...
        mean((rec.g - mean(rec.g)).^2));
    for is = 1:nSeed
        a = rec.a(is:nSeed:end);
        nFrames = numel(a) / p.M;
        nFramesOf(r, is) = nFrames;
        bitCount(r, is) = 4 * p.Nd * numel(a);
        stream = RandStream('mt19937ar', 'Seed', 1000 * r + p.seeds(is));
        for f = 1:nFrames
            af = a((f - 1) * p.M + (1:p.M));
            ab = [af(1); af; af(end)];             % guard blocks take the gain of the adjacent block
            bits = randi(stream, [0 1], 4 * p.Nd, p.B);
            wr = randn(stream, Tmax, 1);
            wi = randn(stream, Tmax, 1);
            w0 = (wr + 1j * wi) / sqrt(2);
            [s0, iq, qq, X] = ofdm_frame(p, bits);
            iq = iq(:, 2:p.M + 1);                 % counted blocks
            qq = qq(:, 2:p.M + 1);
            X = X(:, 2:p.M + 1);
            for filt = 0:1
                fr = prepare_frame(p, ctx, s0, ab, w0, logical(filt));
                for ig = 1:nG
                    [Zs, ~, pOut] = equalized_symbols(p, ctx, fr, p.gammaGrid(ig));
                    alpha = ctx.ref.alpha(ig);
                    pc = ctx.sigma2(filt + 1) * ctx.ref.power(ig);
                    alphaNum(r, is, filt + 1, ig) = alphaNum(r, is, filt + 1, ig) + ...
                        real(sum(alpha * Zs(:) .* conj(X(:))));
                    pSum(r, is, filt + 1, ig) = pSum(r, is, filt + 1, ig) + pOut;
                    eX(r, is, filt + 1, ig) = eX(r, is, filt + 1, ig) + sum(abs(X(:)).^2);
                    eErr(r, is, filt + 1, ig) = eErr(r, is, filt + 1, ig) + sum(abs(Zs(:) - X(:)).^2);
                    kappa = 1 / sqrt(pc);              % transmit scale for unit average power
                    e = count_bit_errors(Zs, fr.Zn, noiseStd / (kappa * alpha), iq, qq);
                    errors(r, filt + 1, ig, :, is) = reshape(errors(r, filt + 1, ig, :, is), 1, nS) + e;
                end
            end
            fprintf('%s run %d frame %d/%d (%.0f s)\n', p.records{r}, is, f, nFrames, toc(tStart));
        end
    end
end
fclose(chanFile);

% A symbol error energy below the rounding level of the arithmetic is recorded as zero.
eErr(eErr < 1e-20 * eX) = 0;

berFile = open_output(fullfile(resultsDir, 'ber.csv'));
fprintf(berFile, 'record,filter,gamma,snr_db,seed,bit_errors,bits\n');
statFile = open_output(fullfile(resultsDir, 'waveform_stats.csv'));
fprintf(statFile, 'record,seed,filter,gamma,sigma2,alpha_ref,alpha_meas,pc_ref,p_meas,e_x,e_err\n');
for r = 1:nR
    for filt = 0:1
        for ig = 1:nG
            for iSnr = 1:nS
                for is = 1:nSeed
                    fprintf(berFile, '%s,%d,%s,%d,%d,%d,%d\n', p.records{r}, filt, ...
                        gamma_text(p.gammaGrid(ig)), p.snrGrid(iSnr), p.seeds(is), ...
                        errors(r, filt + 1, ig, iSnr, is), bitCount(r, is));
                end
            end
        end
    end
    for is = 1:nSeed
        for filt = 0:1
            for ig = 1:nG
                s2 = ctx.sigma2(filt + 1);
                fprintf(statFile, '%s,%d,%d,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
                    p.records{r}, p.seeds(is), filt, gamma_text(p.gammaGrid(ig)), s2, ...
                    ctx.ref.alpha(ig), alphaNum(r, is, filt + 1, ig) / eX(r, is, filt + 1, ig), ...
                    s2 * ctx.ref.power(ig), pSum(r, is, filt + 1, ig) / nFramesOf(r, is), ...
                    eX(r, is, filt + 1, ig), eErr(r, is, filt + 1, ig));
            end
        end
    end
end
fclose(berFile);
fclose(statFile);
end

function fid = open_output(name)
%OPEN_OUTPUT Open a file for writing and stop with an error when that is not possible.
fid = fopen(name, 'w');
if fid < 0
    error('run_simulation:file', 'Cannot write %s.', name);
end
end
