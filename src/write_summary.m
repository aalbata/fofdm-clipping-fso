function write_summary(resultsDir)
%WRITE_SUMMARY Summary of the waveform configurations under both transmitter models.
%   Reads ber.csv, calibration.csv, selected_gamma.csv, and closed_form_reference.csv from
%   resultsDir and writes summary.csv, summary_differences.csv, and peak_gain_sensitivity.csv.
%   Configurations: ofdm and fofdm (conventional OFDM and F-OFDM without a limiter), and
%   each transmit chain at the clipping ratio selected in selected_gamma.csv.
%   summary.csv has one row per transmitter model and configuration, with the required SNR
%   (average-power model) or the required peak SNR (peak-limited model) for each BER target
%   and record rounded to two decimals, the PAPR used by the peak-limited model, and the
%   guard-band level and power beyond half the sampling rate of calibration.csv.
%   summary_differences.csv holds named differences in decibels between entries of
%   summary.csv, between the OFDM entries and the closed form of the same record, and
%   between the closed form of each record and the closed form with unit gain. Each
%   difference is formed from the values rounded to two decimals, so that it equals the
%   difference of the printed entries. The columns run_min_db and run_max_db give the range
%   of the same difference over the individual runs, formed before rounding (NaN for the
%   differences between closed-form values).
p = sim_params();
S = read_ber(resultsDir);
R = required_table(p, S);
[~, ~, cal] = read_calibration(p, resultsDir);
paprDb = chain_papr(p, cal);
closed = readtable(fullfile(resultsDir, 'closed_form_reference.csv'), 'TextType', 'string');
nR = numel(p.records);
nT = numel(p.berTargets);
nSeed = size(R.reqSeed, 5);

cfg = configurations(p, resultsDir);
names = cfg.name;
chainOf = cfg.chain;
gammaOf = cfg.gamma;
nC = numel(names);
models = {'avg', 'peak'};

value = zeros(2, nC, nT, nR);                      % model, configuration, target, record
valueRun = zeros(2, nC, nT, nR, nSeed);
fid = fopen(fullfile(resultsDir, 'summary.csv'), 'w');
fprintf(fid, 'model,configuration,chain,gamma');
for it = 1:nT
    for r = 1:nR
        fprintf(fid, ',req_%s_%s_db', p.records{r}, target_text(p.berTargets(it)));
    end
end
fprintf(fid, ',papr_db,guard_band_db,outside_band_db\n');
for m = 1:2
    for c = 1:nC
        ic = find(strcmp(p.chains, chainOf{c}));
        ig = find(round(10 * p.gammaGrid) == round(10 * gammaOf(c)) | (isinf(p.gammaGrid) & isinf(gammaOf(c))));
        shift = 0;
        if m == 2
            shift = paprDb(ic, ig);                % peak SNR = SNR + PAPR of the transmitted waveform
        end
        row = cal.chain == chainOf{c} & round(10 * cal.gamma) == round(10 * gammaOf(c));
        fprintf(fid, '%s,%s,%s,%s', models{m}, names{c}, chainOf{c}, gamma_text(gammaOf(c)));
        for it = 1:nT
            for r = 1:nR
                value(m, c, it, r) = round(100 * (R.req(r, ic, ig, it) + shift)) / 100;
                valueRun(m, c, it, r, :) = R.reqSeed(r, ic, ig, it, :) + shift;
                fprintf(fid, ',%s', num_text(value(m, c, it, r), 2));
            end
        end
        fprintf(fid, ',%s,%s,%s\n', num_text(paprDb(ic, ig), 2), num_text(cal.guard_band_db(row), 0), ...
            num_text(cal.outside_band_db(row), 2));
    end
end
fclose(fid);

% Each row: name, model and configuration of the minuend, model and configuration of the subtrahend.
idx = @(name) find(strcmp(names, name));
pairs = {
    'fofdm_minus_ofdm_avg',              1, idx('fofdm'),       1, idx('ofdm')
    'clip_minus_ofdm_avg',               1, idx('clip'),        1, idx('ofdm')
    'filter_clip_minus_fofdm_avg',       1, idx('filter_clip'), 1, idx('fofdm')
    'clip_filter_minus_fofdm_avg',       1, idx('clip_filter'), 1, idx('fofdm')
    'rcf_minus_ofdm_avg',                1, idx('rcf'),         1, idx('ofdm')
    'filter_clip_minus_clip_avg',        1, idx('filter_clip'), 1, idx('clip')
    'fofdm_minus_ofdm_peak',             2, idx('fofdm'),       2, idx('ofdm')
    'ofdm_minus_clip_peak',              2, idx('ofdm'),        2, idx('clip')
    'fofdm_minus_filter_clip_peak',      2, idx('fofdm'),       2, idx('filter_clip')
    'fofdm_minus_clip_filter_peak',      2, idx('fofdm'),       2, idx('clip_filter')
    'ofdm_minus_rcf_peak',               2, idx('ofdm'),        2, idx('rcf')
    'ofdm_minus_filter_clip_peak',       2, idx('ofdm'),        2, idx('filter_clip')
    'ofdm_minus_clip_filter_peak',       2, idx('ofdm'),        2, idx('clip_filter')
    'filter_clip_minus_clip_peak',       2, idx('filter_clip'), 2, idx('clip')
    'clip_filter_minus_clip_peak',       2, idx('clip_filter'), 2, idx('clip')
    'rcf_minus_clip_peak',               2, idx('rcf'),         2, idx('clip')
    'clip_filter_minus_filter_clip_peak', 2, idx('clip_filter'), 2, idx('filter_clip')
    };
fid = fopen(fullfile(resultsDir, 'summary_differences.csv'), 'w');
fprintf(fid, 'quantity,ber_target,record,difference_db,run_min_db,run_max_db\n');
for it = 1:nT
    for r = 1:nR
        ref = closed.snr_req_db(closed.ber_target == p.berTargets(it) & closed.channel == p.records{r});
        awgn = closed.snr_req_db(closed.ber_target == p.berTargets(it) & closed.channel == "awgn");
        d = reshape(valueRun(1, idx('ofdm'), it, r, :), 1, []) - ref;
        fprintf(fid, 'ofdm_avg_minus_closed_form,%g,%s,%s,%s,%s\n', p.berTargets(it), p.records{r}, ...
            num_text(value(1, idx('ofdm'), it, r) - round(100 * ref) / 100, 2), ...
            num_text(min(d, [], 'includenan'), 2), num_text(max(d, [], 'includenan'), 2));
        fprintf(fid, 'closed_form_minus_awgn,%g,%s,%s,NaN,NaN\n', p.berTargets(it), p.records{r}, ...
            num_text(round(100 * ref) / 100 - round(100 * awgn) / 100, 2));
        for i = 1:size(pairs, 1)
            d = reshape(valueRun(pairs{i, 2}, pairs{i, 3}, it, r, :) - valueRun(pairs{i, 4}, pairs{i, 5}, it, r, :), 1, []);
            fprintf(fid, '%s,%g,%s,%s,%s,%s\n', pairs{i, 1}, p.berTargets(it), p.records{r}, ...
                num_text(value(pairs{i, 2}, pairs{i, 3}, it, r) - value(pairs{i, 4}, pairs{i, 5}, it, r), 2), ...
                num_text(min(d, [], 'includenan'), 2), num_text(max(d, [], 'includenan'), 2));
        end
    end
end
fclose(fid);

% Peak-power gain of each chain at its selected clipping ratio over its unclipped reference (OFDM
% for clip and rcf, F-OFDM for filter_clip and clip_filter) for the peak power defined at each
% exceedance probability of p.peakProbs, formed from required peak SNR values rounded to two
% decimals as in summary.csv.
qCol = arrayfun(@(q) sprintf('papr_1e%d_db', round(-log10(q))), p.peakProbs, 'UniformOutput', false);
refOf = {'clip', 'clip'; 'filter_clip', 'filter_clip'; 'clip_filter', 'filter_clip'; 'rcf', 'clip'};
fid = fopen(fullfile(resultsDir, 'peak_gain_sensitivity.csv'), 'w');
fprintf(fid, 'chain,gamma,reference,ber_target,record');
fprintf(fid, ',gain_%s', qCol{:});
fprintf(fid, '\n');
iInf = find(isinf(p.gammaGrid));
for i = 1:size(refOf, 1)
    ic = find(strcmp(p.chains, refOf{i, 1}));
    ir = find(strcmp(p.chains, refOf{i, 2}));
    c = find(strcmp(chainOf, refOf{i, 1}) & ~isinf(gammaOf));
    ig = cfg.ig(c);
    for it = 1:nT
        for r = 1:nR
            fprintf(fid, '%s,%s,%s,%g,%s', refOf{i, 1}, gamma_text(gammaOf(c)), refOf{i, 2}, p.berTargets(it), p.records{r});
            for iq = 1:numel(qCol)
                pr = cal.(qCol{iq})(cal.chain == refOf{i, 2} & isinf(cal.gamma));
                pc = cal.(qCol{iq})(cal.chain == refOf{i, 1} & round(10 * cal.gamma) == round(10 * gammaOf(c)));
                g = round(100 * (R.req(r, ir, iInf, it) + pr)) / 100 - round(100 * (R.req(r, ic, ig, it) + pc)) / 100;
                fprintf(fid, ',%s', num_text(g, 2));
            end
            fprintf(fid, '\n');
        end
    end
end
fclose(fid);
end

function t = target_text(target)
%TARGET_TEXT BER target as text for a column name, for example 1e2 for 0.01.
t = sprintf('1e%d', round(-log10(target)));
end
