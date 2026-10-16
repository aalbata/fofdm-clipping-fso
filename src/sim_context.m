function ctx = sim_context(p)
%SIM_CONTEXT Filter responses and limiter constants shared by all frames.
%   Fields: h (subband filter), g (matched receive filter), G (response of the filter
%   cascade on the data subcarriers), phase (phase of the timing advance on the data
%   subcarriers), sigma2 (power at the limiter input without and with the subband filter),
%   ref (Gaussian clipping reference on the clipping-ratio grid).
ctx.h = subband_filter(p);
ctx.g = conj(flipud(ctx.h));
Hf = fft(ctx.h, p.N);
ctx.G = abs(Hf(p.bins)).^2;
ctx.phase = exp(-2j * pi * p.k * p.delta / p.N);
ctx.sigma2 = [limiter_input_power(p, []), limiter_input_power(p, ctx.h)];
ctx.ref = clipping_reference(p.gammaGrid);
end
