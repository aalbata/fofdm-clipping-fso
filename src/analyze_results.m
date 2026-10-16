function analyze_results(dataDir, resultsDir)
%ANALYZE_RESULTS Required SNR, selected clipping ratios, and reference tables.
%   Reads ber.csv and waveform_stats.csv from resultsDir and the power records from dataDir.
%   Writes to resultsDir:
%   required_snr.csv        required SNR (average-power model) and required peak SNR
%                           (peak-limited model) of every configuration
%   selected_gamma.csv      clipping ratio with the lowest required peak SNR
%   clipping_reference.csv  Gaussian clipping reference and simulated in-band SDR
%   validation.csv          simulated BER of OFDM without limiter and the closed form
%   closed_form_reference.csv  required SNR of the closed form without gain fluctuations
%                           and averaged over the gains of each record
%   filter_response.csv     magnitude response of the subband filter
%   Under the peak-limited model the transmit scale is one over the clipping threshold
%   instead of one over the RMS amplitude of the limiter output, so the BER at a peak SNR
%   equals the simulated BER at an SNR lower by the peak-to-average power ratio of the
%   limiter output. The required peak SNR is the required SNR plus that ratio in decibels.
p = sim_params();
S = read_ber(resultsDir);
[req, reqSeed, paprDb] = required_table(p, S);
ref = clipping_reference(p.gammaGrid);
nR = numel(p.records);
nG = numel(p.gammaGrid);
nT = numel(p.berTargets);
finite = ~isinf(p.gammaGrid);

fid = fopen(fullfile(resultsDir, 'required_snr.csv'), 'w');
fprintf(fid, 'record,filter,gamma,ber_target,snr_req_db,snr_req_seed_min_db,snr_req_seed_max_db,papr_db,psnr_req_db\n');
for r = 1:nR
    for iFilt = 1:2
        for ig = 1:nG
            for it = 1:nT
                if finite(ig)
                    papr = paprDb(ig);
                else
                    papr = NaN;                    % the peak-limited model requires a limiter
                end
                fprintf(fid, '%s,%d,%s,%g,%s,%s,%s,%s,%s\n', p.records{r}, iFilt - 1, ...
                    gamma_text(p.gammaGrid(ig)), p.berTargets(it), num_text(req(r, iFilt, ig, it), 4), ...
                    num_text(min(reqSeed(r, iFilt, ig, it, :), [], 'includenan'), 4), ...
                    num_text(max(reqSeed(r, iFilt, ig, it, :), [], 'includenan'), 4), ...
                    num_text(papr, 4), num_text(req(r, iFilt, ig, it) + papr, 4));
            end
        end
    end
end
fclose(fid);

% Clipping ratio that minimizes the required peak SNR at the selection target; of equal
% minima the larger clipping ratio is taken.
isSelection = p.berTargets == p.selectionTarget;
if nnz(isSelection) ~= 1
    error('analyze_results:selectionTarget', 'The selection target must be one of the BER targets.');
end
fid = fopen(fullfile(resultsDir, 'selected_gamma.csv'), 'w');
fprintf(fid, 'record,filter,gamma,psnr_req_db\n');
for r = 1:nR
    for iFilt = 1:2
        psnr = reshape(req(r, iFilt, :, isSelection), 1, nG) + paprDb;
        psnr(~finite) = NaN;
        best = min(psnr);
        ig = find(psnr == best, 1, 'last');
        fprintf(fid, '%s,%d,%s,%s\n', p.records{r}, iFilt - 1, gamma_text(p.gammaGrid(ig)), num_text(best, 4));
    end
end
fclose(fid);

% Gaussian clipping reference and the simulated in-band signal-to-distortion ratio, pooled
% over records and seeds.
W = readtable(fullfile(resultsDir, 'waveform_stats.csv'), 'TextType', 'string');
fid = fopen(fullfile(resultsDir, 'clipping_reference.csv'), 'w');
fprintf(fid, 'gamma,alpha,output_power,distortion,sdr_db,papr_db,clip_probability,sdr_inband_ofdm_db,sdr_inband_fofdm_db\n');
for ig = find(finite)
    gamma = p.gammaGrid(ig);
    rows = round(10 * W.gamma) == round(10 * gamma);
    sdrIn = zeros(1, 2);
    for iFilt = 1:2
        sel = rows & W.filter == iFilt - 1;
        sdrIn(iFilt) = 10 * log10(sum(W.e_x(sel)) / sum(W.e_err(sel)));
    end
    fprintf(fid, '%.1f,%.6f,%.6f,%.6e,%.4f,%.4f,%.6e,%.4f,%.4f\n', gamma, ref.alpha(ig), ref.power(ig), ...
        ref.distortion(ig), 10 * log10(ref.sdr(ig)), paprDb(ig), exp(-gamma^2), sdrIn(1), sdrIn(2));
end
fclose(fid);

% Simulated BER of unfiltered OFDM without limiter against the closed form for 16-QAM,
% averaged over the gains of the record.
fid = fopen(fullfile(resultsDir, 'validation.csv'), 'w');
fprintf(fid, 'record,snr_db,bit_errors,bits,ber_sim,ber_theory\n');
for r = 1:nR
    rec = load_power_record(dataDir, p.records{r});
    for iSnr = 1:numel(p.snrGrid)
        rho = rec.g * (p.N / p.Nd) * 10^(p.snrGrid(iSnr) / 10);
        e = sum(S.errors(r, 1, ~finite, iSnr, :), 'all');
        b = sum(S.bits(r, 1, ~finite, iSnr, :), 'all');
        fprintf(fid, '%s,%d,%d,%d,%.6e,%.6e\n', p.records{r}, p.snrGrid(iSnr), e, b, e / b, mean(qam16_ber(rho)));
    end
end
fclose(fid);

% SNR at which the closed form equals each BER target, without gain fluctuations (unit
% gain) and averaged over the gains of each record.
channels = [{'awgn'}, p.records];
gains = cell(1, nR + 1);
gains{1} = 1;
for r = 1:nR
    rec = load_power_record(dataDir, p.records{r});
    gains{r + 1} = rec.g;
end
fid = fopen(fullfile(resultsDir, 'closed_form_reference.csv'), 'w');
fprintf(fid, 'ber_target,channel,snr_req_db\n');
for it = 1:nT
    for c = 1:numel(channels)
        snr = fzero(@(s) log10(mean(qam16_ber(gains{c} * 10^(s / 10) * p.N / p.Nd))) - log10(p.berTargets(it)), [0 30]);
        fprintf(fid, '%g,%s,%.4f\n', p.berTargets(it), channels{c}, snr);
    end
end
fclose(fid);

% Magnitude response of the subband filter at the subcarrier frequencies k = 0..1600.
Hf = fft(subband_filter(p), p.N);
fid = fopen(fullfile(resultsDir, 'filter_response.csv'), 'w');
fprintf(fid, 'k,magnitude_db\n');
for k = 0:1600
    fprintf(fid, '%d,%s\n', k, num_text(20 * log10(abs(Hf(k + 1))), 4));
end
fclose(fid);
end
