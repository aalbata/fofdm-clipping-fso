function pb = qam16_ber(rho)
%QAM16_BER Bit error rate of Gray-mapped 16-QAM in Gaussian noise.
%   rho is the ratio of the symbol energy to the noise spectral density (linear scale).
x = sqrt(rho / 10);
pb = (3 * erfc(x) + 2 * erfc(3 * x) - erfc(5 * x)) / 8;
end
