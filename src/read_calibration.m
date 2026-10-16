function [txPower, alpha, T] = read_calibration(p, resultsDir)
%READ_CALIBRATION Average transmit power and attenuation factor of every chain and clipping ratio.
%   Reads calibration.csv written by CHARACTERIZE_TRANSMITTERS. txPower and alpha are
%   indexed by chain and clipping ratio in the order of SIM_PARAMS; T is the full table.
T = readtable(fullfile(resultsDir, 'calibration.csv'), 'TextType', 'string');
nC = numel(p.chains);
nG = numel(p.gammaGrid);
txPower = NaN(nC, nG);
alpha = NaN(nC, nG);
for ic = 1:nC
    for ig = 1:nG
        row = T.chain == p.chains{ic} & round(10 * T.gamma) == round(10 * p.gammaGrid(ig));
        if nnz(row) ~= 1
            error('read_calibration:missing', 'calibration.csv has no unique row for %s at gamma %s.', ...
                p.chains{ic}, gamma_text(p.gammaGrid(ig)));
        end
        txPower(ic, ig) = T.tx_power(row);
        alpha(ic, ig) = T.alpha(row);
    end
end
end
