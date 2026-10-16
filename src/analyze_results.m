function analyze_results(dataDir, resultsDir)
%ANALYZE_RESULTS Required SNR and peak SNR, selected clipping ratios, and reference tables.
%   Reads ber.csv and calibration.csv from resultsDir and the power records from dataDir.
%   Writes to resultsDir:
%   required_snr.csv        required SNR (average-power model) and required peak SNR
%                           (peak-limited model) of every record, chain, clipping ratio, and
%                           BER target, with the range over the runs and the required SNR
%                           from the whole-decibel grid and from linear interpolation
%   selected_gamma.csv      for every chain, the clipping ratio with the lowest required peak
%                           SNR at p.selectionTarget for the sunny record, the range of
%                           clipping ratios within p.nearOptimalDb of that minimum, and the
%                           minimizing clipping ratios for the other record and target
%   clipping_reference.csv  Gaussian clipping reference and the measured values of the clip chain
%   validation.csv          simulated BER of OFDM without a limiter and the closed form
%   closed_form_reference.csv  required SNR of the closed form with unit gain and averaged
%                           over the gains of each record
%   filter_response.csv     magnitude response of the subband filter
%   Under the peak-limited model the transmitted waveform is scaled so that its peak power,
%   the instantaneous power exceeded by a fraction p.peakProb of the samples of the
%   reconstructed waveform, equals the peak power limit. Because the noise is additive and
%   independent of the signal, the BER at a peak SNR equals the BER of the average-power
%   model at an SNR lower by the peak-to-average power ratio (PAPR) in calibration.csv; the
%   required peak SNR is the required SNR plus that ratio in decibels.
p = sim_params();
S = read_ber(resultsDir);
R = required_table(p, S);
[~, ~, cal] = read_calibration(p, resultsDir);
nR = numel(p.records);
nC = numel(p.chains);
nG = numel(p.gammaGrid);
nT = numel(p.berTargets);
finite = ~isinf(p.gammaGrid);
paprDb = chain_papr(p, cal);

% Chains that coincide without a limiter must give identical bit errors.
iInf = find(~finite);
same = {'clip', 'rcf'; 'filter_clip', 'clip_filter'};
for i = 1:2
    a = strcmp(p.chains, same{i, 1});
    b = strcmp(p.chains, same{i, 2});
    if ~isequal(S.errors(:, a, iInf, :, :), S.errors(:, b, iInf, :, :))
        error('analyze_results:inf', 'Chains %s and %s differ without a limiter.', same{i, :});
    end
end

fid = fopen(fullfile(resultsDir, 'required_snr.csv'), 'w');
fprintf(fid, ['record,chain,gamma,ber_target,snr_req_db,snr_req_run_min_db,snr_req_run_max_db,' ...
    'snr_req_1db_grid_db,snr_req_linear_db,papr_db,psnr_req_db\n']);
for r = 1:nR
    for ic = 1:nC
        for ig = 1:nG
            for it = 1:nT
                q = R.req(r, ic, ig, it);
                fprintf(fid, '%s,%s,%s,%g,%s,%s,%s,%s,%s,%s,%s\n', p.records{r}, p.chains{ic}, ...
                    gamma_text(p.gammaGrid(ig)), p.berTargets(it), num_text(q, 4), ...
                    num_text(min(R.reqSeed(r, ic, ig, it, :), [], 'includenan'), 4), ...
                    num_text(max(R.reqSeed(r, ic, ig, it, :), [], 'includenan'), 4), ...
                    num_text(R.req1dB(r, ic, ig, it), 4), num_text(R.reqLin(r, ic, ig, it), 4), ...
                    num_text(paprDb(ic, ig), 4), num_text(q + paprDb(ic, ig), 4));
            end
        end
    end
end
fclose(fid);

% Clipping ratio that minimizes the required peak SNR; of equal minima the larger clipping
% ratio is taken. Operation without a limiter is not a candidate.
isSel = find(p.berTargets == p.selectionTarget);
if numel(isSel) ~= 1
    error('analyze_results:selectionTarget', 'The selection target must be one of the BER targets.');
end
iSun = find(strcmp(p.records, 'sun'));
iRain = find(strcmp(p.records, 'rain'));
other = find(p.berTargets ~= p.selectionTarget, 1);
fid = fopen(fullfile(resultsDir, 'selected_gamma.csv'), 'w');
fprintf(fid, ['chain,gamma,psnr_req_db,near_optimal_min,near_optimal_max,near_optimal_count,' ...
    'gamma_rain,gamma_sun_other_target,gamma_rain_other_target\n']);
for ic = 1:nC
    [ig, best, near] = best_gamma(p, R.req(iSun, ic, :, isSel), paprDb(ic, :));
    gNear = p.gammaGrid(near);
    fprintf(fid, '%s,%s,%s,%s,%s,%d,%s,%s,%s\n', p.chains{ic}, gamma_text(p.gammaGrid(ig)), num_text(best, 4), ...
        gamma_text(min(gNear)), gamma_text(max(gNear)), numel(gNear), ...
        gamma_text(p.gammaGrid(best_gamma(p, R.req(iRain, ic, :, isSel), paprDb(ic, :)))), ...
        gamma_text(p.gammaGrid(best_gamma(p, R.req(iSun, ic, :, other), paprDb(ic, :)))), ...
        gamma_text(p.gammaGrid(best_gamma(p, R.req(iRain, ic, :, other), paprDb(ic, :)))));
end
fclose(fid);

% Gaussian clipping reference of the limiter samples and the measured values of the clip
% chain (attenuation factor, in-band SDR, and PAPR of the reconstructed waveform).
ref = clipping_reference(p.gammaGrid);
iClip = find(strcmp(p.chains, 'clip'));
fid = fopen(fullfile(resultsDir, 'clipping_reference.csv'), 'w');
fprintf(fid, ['gamma,alpha,output_power,distortion,sdr_db,papr_db,clip_probability,' ...
    'alpha_clip,inband_sdr_clip_db,papr_1e4_clip_db\n']);
for ig = find(finite)
    gamma = p.gammaGrid(ig);
    row = cal.chain == "clip" & round(10 * cal.gamma) == round(10 * gamma);
    fprintf(fid, '%.1f,%.6f,%.6f,%.6e,%.4f,%.4f,%.6e,%.6f,%.4f,%.4f\n', gamma, ref.alpha(ig), ref.power(ig), ...
        ref.distortion(ig), 10 * log10(ref.sdr(ig)), 10 * log10(ref.papr(ig)), exp(-gamma^2), ...
        cal.alpha(row), cal.inband_sdr_db(row), paprDb(iClip, ig));
end
fclose(fid);

% Simulated BER of unfiltered OFDM without a limiter against the closed form for 16-QAM,
% averaged over the gains of the record.
fid = fopen(fullfile(resultsDir, 'validation.csv'), 'w');
fprintf(fid, 'record,snr_db,bit_errors,bits,ber_sim,ber_theory\n');
for r = 1:nR
    rec = load_power_record(dataDir, p.records{r});
    for iSnr = 1:numel(p.snrGrid)
        rho = rec.g * (p.N / p.Nd) * 10^(p.snrGrid(iSnr) / 10);
        e = sum(S.errors(r, iClip, iInf, iSnr, :), 'all');
        b = sum(S.bits(r, iClip, iInf, iSnr, :), 'all');
        fprintf(fid, '%s,%.2f,%d,%d,%.6e,%.6e\n', p.records{r}, p.snrGrid(iSnr), e, b, e / b, mean(qam16_ber(rho)));
    end
end
fclose(fid);

% SNR at which the closed form equals each BER target, with unit gain and averaged over the
% gains of each record.
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

function [ig, best, near] = best_gamma(p, req, papr)
%BEST_GAMMA Clipping ratio with the lowest required peak SNR and the near-optimal ratios.
psnr = reshape(req, 1, []) + reshape(papr, 1, []);
psnr(isinf(p.gammaGrid)) = NaN;
best = min(psnr);
ig = find(psnr == best, 1, 'last');
near = find(psnr <= best + p.nearOptimalDb);
end
