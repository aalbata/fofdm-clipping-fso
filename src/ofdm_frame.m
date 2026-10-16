function [s0, iq, qq, X] = ofdm_frame(p, bits)
%OFDM_FRAME Map the bits of one frame to 16-QAM and build the cyclic-prefixed block stream.
%   bits is a 4*Nd-by-B array of zeros and ones with one column per block. Within a block,
%   bits 1 and 2 of each group of four select the in-phase level and bits 3 and 4 the
%   quadrature level with the Gray map 00 -> -3, 01 -> -1, 11 -> +1, 10 -> +3.
%   s0 is the column stream of B blocks. iq and qq are the Nd-by-B level indices (0..3)
%   of the in-phase and quadrature components; level = 2*index - 3. X holds the Nd-by-B
%   transmitted symbols with unit average power.
B = size(bits, 2);
b = reshape(bits, 4, p.Nd, B);
pairToIndex = [0 1 3 2];                           % value 2*b1 + b2 of a bit pair -> level index
iq = reshape(pairToIndex(2 * b(1, :, :) + b(2, :, :) + 1), p.Nd, B);
qq = reshape(pairToIndex(2 * b(3, :, :) + b(4, :, :) + 1), p.Nd, B);
X = ((2 * iq - 3) + 1j * (2 * qq - 3)) / sqrt(10);
Xf = zeros(p.N, B);
Xf(p.bins, :) = X;
x = ifft(Xf, [], 1) * sqrt(p.N);                   % unitary inverse transform
blocks = [x(p.N - p.Ncp + 1:p.N, :); x];           % cyclic prefix followed by the block
s0 = blocks(:);
end
