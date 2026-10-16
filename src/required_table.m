function R = required_table(p, S)
%REQUIRED_TABLE Required SNR of every record, chain, clipping ratio, and BER target.
%   S is the structure returned by READ_BER (fields ber and berSeed are used). Fields of R,
%   indexed by record, chain, clipping ratio, and BER target (NaN when the target is not
%   reached):
%   req      SNR in dB at which the pooled BER crosses the target (parabola in log BER)
%   reqSeed  the same for the BER of each run (run as fifth index)
%   req1dB   the same as req from the grid points at whole decibels only
%   reqLin   the same as req with the straight line through log BER (REQUIRED_SNR)
nR = numel(p.records);
nC = numel(p.chains);
nG = numel(p.gammaGrid);
nT = numel(p.berTargets);
nSeed = size(S.berSeed, 5);
whole = abs(p.snrGrid - round(p.snrGrid)) < 1e-9;
R.req = NaN(nR, nC, nG, nT);
R.req1dB = NaN(nR, nC, nG, nT);
R.reqLin = NaN(nR, nC, nG, nT);
R.reqSeed = NaN(nR, nC, nG, nT, nSeed);
for r = 1:nR
    for ic = 1:nC
        for ig = 1:nG
            ber = reshape(S.ber(r, ic, ig, :), 1, []);
            for it = 1:nT
                target = p.berTargets(it);
                R.req(r, ic, ig, it) = required_snr(p.snrGrid, ber, target);
                R.req1dB(r, ic, ig, it) = required_snr(p.snrGrid(whole), ber(whole), target);
                R.reqLin(r, ic, ig, it) = required_snr(p.snrGrid, ber, target, 'linear');
                for is = 1:nSeed
                    R.reqSeed(r, ic, ig, it, is) = required_snr(p.snrGrid, ...
                        reshape(S.berSeed(r, ic, ig, :, is), 1, []), target);
                end
            end
        end
    end
end
end
