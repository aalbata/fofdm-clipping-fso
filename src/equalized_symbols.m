function [Zs, A, pOut] = equalized_symbols(p, ctx, fr, gamma)
%EQUALIZED_SYMBOLS Noise-free equalized symbols of the counted blocks for one clipping ratio.
%   fr is the frame prepared by PREPARE_FRAME and gamma the clipping ratio (Inf for no
%   limiter). The limiter acts on the stream interpolated by the factor p.os; the receiver
%   observes the limiter output within the sampling bandwidth (ideal lowpass filter and
%   decimation). Zs is Nd-by-M: the FFT outputs divided by the known response, which
%   includes the attenuation factor of the limiter. A is the clipping threshold and pOut
%   the mean power of the limiter output over the counted blocks. The transmit scale
%   factor cancels in Zs and is therefore omitted.
if isinf(gamma)
    A = Inf;
    alpha = 1;
    c = fr.u;
    pOut = mean(abs(fr.uo(fr.powerRange)).^2);
else
    A = gamma * fr.sigma;
    ref = clipping_reference(gamma);
    alpha = ref.alpha;
    co = envelope_limiter(fr.uo, A);
    pOut = mean(abs(co(fr.powerRange)).^2);
    c = resample_stream(co, p.nPad);
    c = c(1:numel(fr.u));
end
y = fr.an .* c;
if fr.filt
    z = filter_stream(y, ctx.g);
else
    z = y;
end
Zs = receive_frame(p, z, fr.Dc) ./ (alpha * fr.eqz);
end
