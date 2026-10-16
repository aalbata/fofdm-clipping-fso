function e = count_bit_errors(Zs, Zn, scale, iq, qq)
%COUNT_BIT_ERRORS Bit errors of hard 16-QAM decisions for a set of noise scales.
%   Decisions are taken on Zs + scale(i)*Zn, where Zs holds the equalized noise-free
%   symbols and Zn the equalized noise for unit noise scale. iq and qq are the transmitted
%   level indices (0..3) of the in-phase and quadrature components. e(i) is the number of
%   bit errors for scale(i) >= 0.
%   Each component carries two bits. With the component scaled so that the levels are
%   -1.5, -0.5, 0.5, 1.5, the first bit is one for positive values and the second bit is
%   one for magnitudes below one. For a component v0 + s*n0 the set of noise scales s at
%   which a bit is in error is bounded by the scales at which the component crosses 0, -1,
%   or +1, so the errors at all noise scales follow from counts of these crossing scales.
s0 = sqrt(10) / 2;
v0 = [real(Zs(:)); imag(Zs(:))] * s0;
n0 = [real(Zn(:)); imag(Zn(:))] * s0;
level = [iq(:); qq(:)];
sentPositive = level >= 2;
sentInner = level == 1 | level == 2;
scale = scale(:);

% Components without noise: errors do not depend on the noise scale.
zero = n0 == 0;
base = nnz(xor(v0(zero) > 0, sentPositive(zero))) + nnz(xor(abs(v0(zero)) < 1, sentInner(zero)));
v = v0(~zero);
n = n0(~zero);
sp = sentPositive(~zero);
si = sentInner(~zero);

% First bit: for n > 0 the component is positive for s > t, for n < 0 for s < t.
t = -v ./ n;
up = n > 0;
e1 = nnz(up & sp) - count_less(t(up & sp), scale) ...          % s <= t
    + count_less(t(up & ~sp), scale) ...                      % s > t
    + count_less_equal(t(~up & sp), scale) ...                % s >= t
    + nnz(~up & ~sp) - count_less_equal(t(~up & ~sp), scale); % s < t

% Second bit: the magnitude is below one for l < s < u.
t1 = (-1 - v) ./ n;
t2 = (1 - v) ./ n;
l = min(t1, t2);
u = max(t1, t2);
e2 = nnz(si) - count_less(l(si), scale) + count_less_equal(u(si), scale) ...   % s <= l or s >= u
    + count_less(l(~si), scale) - count_less_equal(u(~si), scale);             % l < s < u
e = reshape(base + e1 + e2, 1, []);
end

function c = count_less(a, s)
%COUNT_LESS Number of elements of a below each element of s.
[ss, order] = sort(s);
h = histcounts(a, [-Inf; ss; Inf]);
cs = cumsum(h(:));
c = zeros(numel(s), 1);
c(order) = cs(1:numel(s));
end

function c = count_less_equal(a, s)
%COUNT_LESS_EQUAL Number of elements of a at or below each element of s.
c = numel(a) - count_less(-a, -s);
end
