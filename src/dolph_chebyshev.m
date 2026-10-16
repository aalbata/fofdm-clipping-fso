function w = dolph_chebyshev(L, r)
%DOLPH_CHEBYSHEV Dolph-Chebyshev window of length L with sidelobes r dB below the main lobe.
%   The discrete Fourier transform of the window, taken about its center, is the
%   Chebyshev polynomial of degree L-1 evaluated at beta*cos(pi*k/L), k = 0..L-1, where
%   beta sets the ratio of the main lobe to the sidelobes to 10^(r/20). The window is the
%   corresponding cosine sum, normalized to a maximum of one.
degree = L - 1;
beta = cosh(acosh(10^(r / 20)) / degree);
k = 0:L - 1;
x = beta * cos(pi * k / L);
c = zeros(1, L);
inside = abs(x) <= 1;
c(inside) = cos(degree * acos(x(inside)));
c(~inside) = sign(x(~inside)).^degree .* cosh(degree * acosh(abs(x(~inside))));
n = (0:L - 1).' - (L - 1) / 2;                     % sample positions relative to the center
w = cos(2 * pi * n * k / L) * c.';
w = w / max(w);
end
