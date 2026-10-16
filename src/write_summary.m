function write_summary(resultsDir)
%WRITE_SUMMARY Summary table of the four waveform configurations under both transmitter models.
%   Reads ber.csv, selected_gamma.csv, transmit_stats.csv, and closed_form_reference.csv
%   from resultsDir and writes summary.csv and summary_differences.csv.
%   summary.csv has one row per transmitter model and configuration, with the required SNR
%   (average-power model) or the required peak SNR (peak-limited model) for each BER target
%   and record, and the guard-band level of the transmitted stream. Under the average-power
%   model OFDM and F-OFDM run without limiter; under the peak-limited model they run with
%   the clipping ratio p.gamma0. Clipped OFDM and joint clipping-assisted F-OFDM use the
%   clipping ratio selected on the sunny record under both models.
%   summary_differences.csv holds differences in decibels between entries of summary.csv,
%   between the OFDM entries and the closed form of the same record, and between the closed
%   form of each record and the closed form with unit gain. Each difference is formed from
%   values rounded to two decimals. The columns seed_min_db and seed_max_db give the range of
%   the same difference over the individual runs, formed from the per-run values before
%   rounding (NaN for the differences between closed-form values).
p = sim_params();
st = figure_style();
S = read_ber(resultsDir);
[req, reqSeed, paprDb] = required_table(p, S);
sel = readtable(fullfile(resultsDir, 'selected_gamma.csv'), 'TextType', 'string');
stats = readtable(fullfile(resultsDir, 'transmit_stats.csv'), 'TextType', 'string');
closed = readtable(fullfile(resultsDir, 'closed_form_reference.csv'), 'TextType', 'string');
gammaSel = sel.gamma(sel.record == "sun" & sel.filter == 0);
if gammaSel ~= sel.gamma(sel.record == "sun" & sel.filter == 1)
    error('write_summary:selection', 'The selected clipping ratio differs between the filter settings.');
end
nR = numel(p.records);
nT = numel(p.berTargets);
nSeed = numel(p.seeds);

models = {'avg', 'peak'};
filt = [0 1 0 1];
gammaOf = {[Inf, Inf, gammaSel, gammaSel], [p.gamma0, p.gamma0, gammaSel, gammaSel]};
guardName = {{'ofdm', 'fofdm', 'clipped_ofdm', 'joint'}, {'ofdm_gamma0', 'fofdm_gamma0', 'clipped_ofdm', 'joint'}};

value = zeros(2, 4, nT, nR);                       % model, configuration, target, record
valueSeed = zeros(2, 4, nT, nR, nSeed);
fid = fopen(fullfile(resultsDir, 'summary.csv'), 'w');
fprintf(fid, 'model,configuration,filter,gamma');
for it = 1:nT
    for r = 1:nR
        fprintf(fid, ',req_%s_%s_db', p.records{r}, target_text(p.berTargets(it)));
    end
end
fprintf(fid, ',guard_band_db\n');
for m = 1:2
    for c = 1:4
        gamma = gammaOf{m}(c);
        ig = find(round(10 * p.gammaGrid) == round(10 * gamma));
        shift = 0;
        if m == 2
            shift = paprDb(ig);                    % peak SNR = SNR + PAPR of the limiter output
        end
        fprintf(fid, '%s,%s,%d,%s', models{m}, st.names{c}, filt(c), gamma_text(gamma));
        for it = 1:nT
            for r = 1:nR
                value(m, c, it, r) = round(100 * (req(r, filt(c) + 1, ig, it) + shift)) / 100;
                valueSeed(m, c, it, r, :) = reqSeed(r, filt(c) + 1, ig, it, :) + shift;
                fprintf(fid, ',%s', num_text(value(m, c, it, r), 2));
            end
        end
        fprintf(fid, ',%s\n', num_text(stats.guard_band_db(stats.configuration == guardName{m}{c}), 0));
    end
end
fclose(fid);

% Each row: name, model and configuration of the minuend, model and configuration of the subtrahend.
pairs = {
    'fofdm_minus_ofdm_avg',           1, 2, 1, 1
    'clipped_ofdm_minus_ofdm_avg',    1, 3, 1, 1
    'joint_minus_fofdm_avg',          1, 4, 1, 2
    'joint_minus_clipped_ofdm_avg',   1, 4, 1, 3
    'ofdm_minus_clipped_ofdm_peak',   2, 1, 2, 3
    'fofdm_minus_joint_peak',         2, 2, 2, 4
    'joint_minus_clipped_ofdm_peak',  2, 4, 2, 3
    };
fid = fopen(fullfile(resultsDir, 'summary_differences.csv'), 'w');
fprintf(fid, 'quantity,ber_target,record,difference_db,seed_min_db,seed_max_db\n');
for it = 1:nT
    for r = 1:nR
        ref = closed.snr_req_db(closed.ber_target == p.berTargets(it) & closed.channel == p.records{r});
        awgn = closed.snr_req_db(closed.ber_target == p.berTargets(it) & closed.channel == "awgn");
        d = reshape(valueSeed(1, 1, it, r, :), 1, []) - ref;
        fprintf(fid, 'ofdm_avg_minus_closed_form,%g,%s,%s,%s,%s\n', p.berTargets(it), p.records{r}, ...
            num_text(value(1, 1, it, r) - round(100 * ref) / 100, 2), num_text(min(d, [], 'includenan'), 2), ...
            num_text(max(d, [], 'includenan'), 2));
        fprintf(fid, 'closed_form_minus_awgn,%g,%s,%s,NaN,NaN\n', p.berTargets(it), p.records{r}, ...
            num_text(round(100 * ref) / 100 - round(100 * awgn) / 100, 2));
        for i = 1:size(pairs, 1)
            d = reshape(valueSeed(pairs{i, 2}, pairs{i, 3}, it, r, :) - valueSeed(pairs{i, 4}, pairs{i, 5}, it, r, :), 1, []);
            fprintf(fid, '%s,%g,%s,%s,%s,%s\n', pairs{i, 1}, p.berTargets(it), p.records{r}, ...
                num_text(value(pairs{i, 2}, pairs{i, 3}, it, r) - value(pairs{i, 4}, pairs{i, 5}, it, r), 2), ...
                num_text(min(d, [], 'includenan'), 2), num_text(max(d, [], 'includenan'), 2));
        end
    end
end
fclose(fid);
end

function t = target_text(target)
%TARGET_TEXT BER target as text for a column name, for example 1e2 for 0.01.
t = sprintf('1e%d', round(-log10(target)));
end
