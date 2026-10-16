function [xs, xn, rate] = transmit_stream(p, ctx, fr0, frF, chain, gamma, X)
%TRANSMIT_STREAM Transmitted waveform of one frame for a transmit chain and clipping ratio.
%   fr0 and frF are the frames prepared by PREPARE_FRAME without and with the subband
%   filter, chain is one of p.chains, gamma the clipping ratio (Inf for no limiter), and X
%   the Nd-by-B data symbols of the frame. The clipping threshold is gamma times the RMS
%   amplitude of the limiter input.
%   xn is the transmitted waveform at rate times the sampling rate, as the transmitter
%   produces it: the limiter output at p.os times the sampling rate for clip and
%   filter_clip, and a stream at the sampling rate (rate = 1) for clip_filter, rcf, and
%   operation without a limiter. Its first sample is aligned with the first sample of the
%   prepared frame of the same filter setting, and it is zero-padded to rate * p.nPad samples.
%   xs is the waveform observed within the sampling bandwidth, sampled at the sampling rate,
%   with the length of the prepared frame (fr0.u for clip and rcf, frF.u otherwise).
ic = find(strcmp(p.chains, chain));
A = gamma * sqrt(ctx.chainSigma2(ic));
filtered = p.chainFiltered(ic);
if filtered
    fr = frF;
else
    fr = fr0;
end
T = numel(fr.u);
if isinf(gamma)
    xs = fr.u;
    xn = [fr.u; zeros(p.nPad - T, 1)];
    rate = 1;
    return
end
switch chain
    case {'clip', 'filter_clip'}
        xn = envelope_limiter(fr.uo, A);
        rate = p.os;
        c = resample_stream(xn, p.nPad);
        xs = c(1:T);
    case 'clip_filter'
        c = resample_stream(envelope_limiter(fr0.uo, A), p.nPad);
        xs = filter_stream(c(1:numel(fr0.u)), ctx.h);
        xn = [xs; zeros(p.nPad - T, 1)];
        rate = 1;
    case 'rcf'
        xs = rcf_stream(p, X, A);
        xn = [xs; zeros(p.nPad - T, 1)];
        rate = 1;
    otherwise
        error('transmit_stream:chain', 'Unknown transmit chain "%s".', chain);
end
end
