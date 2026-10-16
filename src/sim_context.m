function ctx = sim_context(p)
%SIM_CONTEXT Filter responses and limiter input powers shared by all frames.
%   Fields: h (subband filter), g (matched receive filter), G (response of the filter
%   cascade on the data subcarriers), phase (phase of the timing advance on the data
%   subcarriers), sigma2 (power of the OFDM stream without and with the subband filter),
%   chainSigma2 (power at the limiter input of each transmit chain).
ctx.h = subband_filter(p);
ctx.g = conj(flipud(ctx.h));
Hf = fft(ctx.h, p.N);
ctx.G = abs(Hf(p.bins)).^2;
ctx.phase = exp(-2j * pi * p.k * p.delta / p.N);
ctx.sigma2 = [limiter_input_power(p, []), limiter_input_power(p, ctx.h)];
ctx.chainSigma2 = ctx.sigma2(1) * ones(1, numel(p.chains));  % limiter input power of each chain
ctx.chainSigma2(strcmp(p.chains, 'filter_clip')) = ctx.sigma2(2);
end
