function t = gamma_text(gamma)
%GAMMA_TEXT Clipping ratio as text with one decimal, or Inf for no limiter.
if isinf(gamma)
    t = 'Inf';
else
    t = sprintf('%.1f', gamma);
end
end
