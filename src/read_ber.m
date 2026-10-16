function S = read_ber(resultsDir)
%READ_BER Bit-error counts of ber.csv as arrays.
%   Fields: errors and bits, indexed by record, chain, clipping ratio, SNR, and run in the
%   order of SIM_PARAMS; ber, the BER pooled over the runs; berSeed, the BER of each run.
p = sim_params();
T = readtable(fullfile(resultsDir, 'ber.csv'), 'TextType', 'string');
nSeed = numel(p.seeds);
[~, r] = ismember(T.record, string(p.records));
[~, ic] = ismember(T.chain, string(p.chains));
[~, ig] = ismember(round(10 * T.gamma), round(10 * p.gammaGrid));
[~, iSnr] = ismember(round(100 * T.snr_db), round(100 * p.snrGrid));
if any([r; ic; ig; iSnr] == 0)
    error('read_ber:grid', 'ber.csv does not match the grids of sim_params.');
end
sz = [numel(p.records), numel(p.chains), numel(p.gammaGrid), numel(p.snrGrid), nSeed];
S.errors = zeros(sz);
S.bits = zeros(sz);
for is = 1:nSeed
    idx = sub2ind(sz, r, ic, ig, iSnr, is * ones(size(r)));
    S.errors(idx) = T.(sprintf('errors_run%d', p.seeds(is)));
    S.bits(idx) = T.bits_per_run;
end
S.ber = sum(S.errors, 5) ./ sum(S.bits, 5);
S.berSeed = S.errors ./ S.bits;
end
