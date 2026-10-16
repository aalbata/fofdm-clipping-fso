function c = envelope_limiter(u, A)
%ENVELOPE_LIMITER Limit the magnitude of complex samples to A and keep their phase.
c = u;
mag = abs(u);
over = mag > A;
c(over) = A * u(over) ./ mag(over);
end
