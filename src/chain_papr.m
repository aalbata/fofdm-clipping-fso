function paprDb = chain_papr(p, cal)
%CHAIN_PAPR Peak-to-average power ratio of every chain and clipping ratio from calibration.csv.
%   cal is the table returned by READ_CALIBRATION. The ratio is the instantaneous power
%   exceeded by a fraction p.peakProb of the samples of the reconstructed waveform relative
%   to its average power, in decibels, indexed by chain and clipping ratio.
if p.peakProb ~= 1e-4
    error('chain_papr:peakProb', 'calibration.csv holds the peak power at an exceedance of 1e-4.');
end
nC = numel(p.chains);
nG = numel(p.gammaGrid);
paprDb = NaN(nC, nG);
for ic = 1:nC
    for ig = 1:nG
        row = cal.chain == p.chains{ic} & round(10 * cal.gamma) == round(10 * p.gammaGrid(ig));
        paprDb(ic, ig) = cal.papr_1e4_db(row);
    end
end
end
