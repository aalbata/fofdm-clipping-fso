function S = read_ber(resultsDir)
%READ_BER Bit-error counts of ber.csv as arrays.
%   Fields: errors and bits, indexed by record, filter (1 = off, 2 = on), clipping ratio,
%   SNR, and seed in the order of the grids of SIM_PARAMS; ber, the BER pooled over the
%   seeds; berSeed, the BER of each seed.
p = sim_params();
T = readtable(fullfile(resultsDir, 'ber.csv'), 'TextType', 'string');
[~, r] = ismember(T.record, string(p.records));
[~, ig] = ismember(round(10 * T.gamma), round(10 * p.gammaGrid));
[~, iSnr] = ismember(T.snr_db, p.snrGrid);
[~, is] = ismember(T.seed, p.seeds);
sz = [numel(p.records), 2, numel(p.gammaGrid), numel(p.snrGrid), numel(p.seeds)];
idx = sub2ind(sz, r, T.filter + 1, ig, iSnr, is);
S.errors = zeros(sz);
S.errors(idx) = T.bit_errors;
S.bits = zeros(sz);
S.bits(idx) = T.bits;
S.ber = sum(S.errors, 5) ./ sum(S.bits, 5);
S.berSeed = S.errors ./ S.bits;
end
