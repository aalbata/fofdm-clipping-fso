function t = num_text(x, decimals)
%NUM_TEXT Number as text with a fixed number of decimals; NaN and infinite values as NaN, Inf, -Inf.
%   A value that rounds to zero is written without a minus sign.
if isnan(x)
    t = 'NaN';
elseif isinf(x)
    if x > 0
        t = 'Inf';
    else
        t = '-Inf';
    end
else
    scale = 10^decimals;
    t = sprintf('%.*f', decimals, round(x * scale) / scale + 0);
end
end
