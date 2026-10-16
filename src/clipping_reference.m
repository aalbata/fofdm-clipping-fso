function ref = clipping_reference(gamma)
%CLIPPING_REFERENCE Gaussian reference of the envelope limiter for clipping ratio gamma.
%   gamma is the ratio of the clipping threshold to the RMS amplitude of a zero-mean
%   circularly symmetric complex Gaussian input. All powers are relative to the input power.
%   Fields: alpha (attenuation factor), power (output power), distortion (distortion
%   power), sdr (signal-to-distortion ratio), papr (peak-to-average power ratio of the
%   output). gamma = Inf describes operation without a limiter.
ref.alpha = 1 - exp(-gamma.^2) + sqrt(pi) / 2 * gamma .* erfc(gamma);
ref.alpha(isinf(gamma)) = 1;
ref.power = 1 - exp(-gamma.^2);
ref.distortion = ref.power - ref.alpha.^2;
ref.sdr = ref.alpha.^2 ./ ref.distortion;
ref.papr = gamma.^2 ./ ref.power;
end
