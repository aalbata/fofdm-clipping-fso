function level = exceedance_level(counts, edges, total, q)
%EXCEEDANCE_LEVEL Level exceeded by a fraction q of the samples, from a histogram of levels.
%   counts(j) is the number of samples with level in [edges(j), edges(j+1)), total the
%   number of samples including those below edges(1), and q a vector of exceedance
%   fractions. The level is interpolated linearly within the bin in which the number of
%   samples above it reaches q*total. NaN is returned when fewer than q*total samples lie
%   above edges(1) or when no sample is counted at all.
counts = counts(:);
edges = edges(:);
above = flipud(cumsum(flipud(counts)));            % samples at or above edges(j)
above = [above; 0];
level = NaN(size(q));
for i = 1:numel(q)
    target = q(i) * total;
    j = find(above >= target, 1, 'last');          % edges(j) has at least target samples above
    if isempty(j) || j > numel(counts) || target <= 0
        continue
    end
    % Within bin j the count above decreases linearly from above(j) to above(j + 1).
    frac = (above(j) - target) / max(above(j) - above(j + 1), eps);
    level(i) = edges(j) + frac * (edges(j + 1) - edges(j));
end
end
