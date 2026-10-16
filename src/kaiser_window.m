function w = kaiser_window(n, beta)
%KAISER_WINDOW Symmetric Kaiser window of n samples with shape parameter beta.
x = (2 * (0:n - 1).' - (n - 1)) / (n - 1);
w = besseli(0, beta * sqrt(1 - x.^2)) / besseli(0, beta);
end
