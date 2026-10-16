function cfg = configurations(p, resultsDir)
%CONFIGURATIONS The six configurations of the summary, the tables, and the figures.
%   OFDM and F-OFDM without a limiter, and each transmit chain at the clipping ratio selected
%   in selected_gamma.csv. Fields: name, chain (transmit chain of SIM_PARAMS), gamma, ic and
%   ig (indices of the chain and of the clipping ratio in SIM_PARAMS).
sel = readtable(fullfile(resultsDir, 'selected_gamma.csv'), 'TextType', 'string');
cfg.name = [{'ofdm', 'fofdm'}, p.chains];
cfg.chain = [{'clip', 'filter_clip'}, p.chains];
cfg.gamma = Inf(1, numel(cfg.name));
for ic = 1:numel(p.chains)
    cfg.gamma(2 + ic) = sel.gamma(sel.chain == p.chains{ic});
end
cfg.ic = zeros(1, numel(cfg.name));
cfg.ig = zeros(1, numel(cfg.name));
for c = 1:numel(cfg.name)
    cfg.ic(c) = find(strcmp(p.chains, cfg.chain{c}));
    cfg.ig(c) = find(round(10 * p.gammaGrid) == round(10 * cfg.gamma(c)));
end
end
