function s2 = limiter_input_power(p, h)
%LIMITER_INPUT_POWER Time-averaged power of the block stream at the limiter input.
%   h is the impulse response of the subband filter, or [] when the stream is not filtered.
%   The data symbols have unit power and are independent across subcarriers and blocks.
if isempty(h)
    s2 = p.Nd / p.N;
    return
end
m = (-(p.L - 1):(p.L - 1)).';
rh = conv(h, flipud(h));                        % autocorrelation of h at lags m
D = sum(cos(2 * pi * m * p.k.' / p.N), 2);      % sum of exp(j 2 pi k m / N) over the data subcarriers
s2 = sum((p.Nb - abs(m)) .* rh .* D) / (p.Nb * p.N);
end
