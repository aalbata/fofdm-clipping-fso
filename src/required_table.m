function [req, reqSeed, paprDb] = required_table(p, S)
%REQUIRED_TABLE Required SNR of every configuration and the PAPR of its limiter output.
%   S is the structure returned by READ_BER. req is indexed by record, filter (1 = off,
%   2 = on), clipping ratio, and BER target and holds the SNR in dB at which the pooled BER
%   crosses the target (NaN when it does not). reqSeed has the seed as fifth index and holds
%   the same quantity for the BER of each seed. paprDb is the peak-to-average power ratio
%   of the limiter output in dB for each clipping ratio (Inf without limiter). Under the
%   peak-limited model the required peak SNR in dB is req + paprDb.
ref = clipping_reference(p.gammaGrid);
paprDb = 10 * log10(ref.papr);
nR = numel(p.records);
nG = numel(p.gammaGrid);
nT = numel(p.berTargets);
nSeed = numel(p.seeds);
req = NaN(nR, 2, nG, nT);
reqSeed = NaN(nR, 2, nG, nT, nSeed);
for r = 1:nR
    for iFilt = 1:2
        for ig = 1:nG
            for it = 1:nT
                req(r, iFilt, ig, it) = required_snr(p.snrGrid, squeeze(S.ber(r, iFilt, ig, :)), p.berTargets(it));
                for is = 1:nSeed
                    reqSeed(r, iFilt, ig, it, is) = required_snr(p.snrGrid, ...
                        squeeze(S.berSeed(r, iFilt, ig, :, is)), p.berTargets(it));
                end
            end
        end
    end
end
end
