function e = count_bit_errors(Zs, Zn, scale, iq, qq)
%COUNT_BIT_ERRORS Bit errors of hard 16-QAM decisions for a set of noise scales.
%   Decisions are taken on Zs + scale(i)*Zn, where Zs holds the equalized noise-free
%   symbols and Zn the equalized noise for unit noise scale. iq and qq are the transmitted
%   level indices (0..3) of the in-phase and quadrature components. e(i) is the number of
%   bit errors for scale(i).
%   Each component carries two bits. With the component scaled so that the levels are
%   -1.5, -0.5, 0.5, 1.5, the first bit is one for positive values and the second bit is
%   one for magnitudes below one.
s = sqrt(10) / 2;
v0 = [real(Zs(:)); imag(Zs(:))] * s;
n0 = [real(Zn(:)); imag(Zn(:))] * s;
level = [iq(:); qq(:)];
sentPositive = level >= 2;
sentInner = level == 1 | level == 2;
e = zeros(1, numel(scale));
for i = 1:numel(scale)
    v = v0 + scale(i) * n0;
    e(i) = nnz(xor(v > 0, sentPositive)) + nnz(xor(abs(v) < 1, sentInner));
end
end
