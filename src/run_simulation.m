function run_simulation(dataDir, resultsDir)
%RUN_SIMULATION Bit-error simulation of all transmit chains over both power records.
%   Reads the power records from dataDir and calibration.csv from resultsDir, and writes
%   ber.csv and channel_stats.csv to resultsDir. The transmit power of every chain and
%   clipping ratio is normalized to one with the average power in calibration.csv, so the
%   noise variance for an SNR of x dB is 10^(-x/10). The receiver divides by the attenuation
%   factor in calibration.csv. Every configuration of a frame uses the same bits, channel
%   gains, and noise realization. Run i of a record uses the record samples i, i + nSeed,
%   i + 2*nSeed, ..., one sample per counted block, so that the runs together use every
%   sample once.
%   ber.csv has one row per record, chain, clipping ratio, and SNR, with the number of
%   bits of each run and the bit errors of each run.
p = sim_params();
ctx = sim_context(p);
nR = numel(p.records);
nC = numel(p.chains);
nG = numel(p.gammaGrid);
nS = numel(p.snrGrid);
nSeed = numel(p.seeds);
noiseStd = 10.^(-p.snrGrid / 20);                  % noise standard deviation for unit transmit power
Tmax = p.B * p.Nb + p.L - 1;
[txPower, alpha] = read_calibration(p, resultsDir);

errors = zeros(nR, nC, nG, nS, nSeed);             % record, chain, gamma, SNR, run
bitCount = zeros(nR, nSeed);
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
            fr = {prepare_frame(p, ctx, s0, ab, w0, false), prepare_frame(p, ctx, s0, ab, w0, true)};
            for ic = 1:nC
                frc = fr{p.chainFiltered(ic) + 1};
                for ig = 1:nG
                    xs = transmit_stream(p, ctx, fr{1}, fr{2}, p.chains{ic}, p.gammaGrid(ig), X);
                    Zs = equalized_symbols(p, ctx, frc, xs) / alpha(ic, ig);
                    scale = noiseStd * sqrt(txPower(ic, ig)) / alpha(ic, ig);
                    e = count_bit_errors(Zs, frc.Zn, scale, iq, qq);
                    errors(r, ic, ig, :, is) = reshape(errors(r, ic, ig, :, is), 1, nS) + e;
                end
            end
            fprintf('%s run %d frame %d/%d (%.0f s)\n', p.records{r}, is, f, nFrames, toc(tStart));
        end
    end
end
fclose(chanFile);

berFile = open_output(fullfile(resultsDir, 'ber.csv'));
fprintf(berFile, 'record,chain,gamma,snr_db,bits_per_run');
fprintf(berFile, ',errors_run%d', p.seeds);
fprintf(berFile, '\n');
for r = 1:nR
    if any(bitCount(r, :) ~= bitCount(r, 1))
        error('run_simulation:bits', 'The runs of record %s differ in length.', p.records{r});
    end
    for ic = 1:nC
        for ig = 1:nG
            for iSnr = 1:nS
                fprintf(berFile, '%s,%s,%s,%.2f,%d', p.records{r}, p.chains{ic}, gamma_text(p.gammaGrid(ig)), ...
                    p.snrGrid(iSnr), bitCount(r, 1));
                fprintf(berFile, ',%d', errors(r, ic, ig, iSnr, :));
                fprintf(berFile, '\n');
            end
        end
    end
end
fclose(berFile);
end

function fid = open_output(name)
%OPEN_OUTPUT Open a file for writing and stop with an error when that is not possible.
fid = fopen(name, 'w');
if fid < 0
    error('run_simulation:file', 'Cannot write %s.', name);
end
end
