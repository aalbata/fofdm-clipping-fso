function rec = load_power_record(dataDir, name)
%LOAD_POWER_RECORD Received-power record with its unit-mean power gain and amplitude gain.
%   name is 'sun' or 'rain'. Fields: time (s), power (dBm), g (power gain with unit mean),
%   a (amplitude gain, the square root of g).
T = readmatrix(fullfile(dataDir, ['power_' name '.csv']), 'NumHeaderLines', 1);
rec.time = T(:, 1);
rec.power = T(:, 2);
lin = 10.^(rec.power / 10);
rec.g = lin / mean(lin);
rec.a = sqrt(rec.g);
end
