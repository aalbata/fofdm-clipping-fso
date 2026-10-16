function fr = prepare_frame(p, ctx, s0, ab, w0, filt)
%PREPARE_FRAME Limiter input, channel gains, and equalized noise of one frame.
%   s0 is the block stream, ab the amplitude gains of the B blocks, w0 unit-variance
%   complex noise, and filt selects the subband filter (true) or no filter (false).
%   Fields of fr: u (limiter input at the sampling rate), uo (u zero-padded to p.nPad
%   samples and interpolated by the oversampling factor p.os of the limiter), an (gain of
%   each sample of u), D and Dc (delays of the transmit filter and of the filter cascade),
%   eqz (known response on the data subcarriers of the counted blocks, without the
%   transmit scale and the limiter attenuation), Zn (noise at the FFT output divided by
%   eqz), filt (true when the subband filter is used). The transmit filter delay D is
%   L/2 samples, half a sample more than the group delay of the filter.
if filt
    fr.u = filter_stream(s0, ctx.h);
    fr.D = p.L / 2;
    fr.Dc = p.L - 1;
    Gk = ctx.G;
    zn = filter_stream(w0(1:numel(fr.u)), ctx.g);
else
    fr.u = s0;
    fr.D = 0;
    fr.Dc = 0;
    Gk = ones(p.Nd, 1);
    zn = w0(1:numel(s0));
end
fr.filt = filt;
T = numel(fr.u);
fr.uo = resample_stream([fr.u; zeros(p.nPad - T, 1)], p.os * p.nPad);
n = (0:T - 1).';
blockOf = min(max(floor((n - fr.D) / p.Nb), 0), p.B - 1) + 1;
fr.an = ab(blockOf);
fr.eqz = (Gk .* ctx.phase) .* ab(2:p.M + 1).';
fr.Zn = receive_frame(p, zn, fr.Dc) ./ fr.eqz;
end
