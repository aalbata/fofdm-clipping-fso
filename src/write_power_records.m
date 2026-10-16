function write_power_records(sourceDir, dataDir)
%WRITE_POWER_RECORDS Convert the measurement files of the dataset to CSV files in dataDir.
%   Each output file has the columns time_s (time since the first sample, in seconds) and
%   power_dbm (received optical power, at the 0.01 dBm resolution of the measurement).
files = {'2019-04-10_powerMeas_sun.mat', 'sun'; '2019-04-15_powerMeas_rain.mat', 'rain'};
for i = 1:size(files, 1)
    S = load(fullfile(sourceDir, files{i, 1}), 'power', 'samplingTime');
    power = round(100 * S.power(:)) / 100;
    time = (0:numel(power) - 1).' * S.samplingTime;
    fid = fopen(fullfile(dataDir, ['power_' files{i, 2} '.csv']), 'w');
    fprintf(fid, 'time_s,power_dbm\n');
    fprintf(fid, '%d,%.2f\n', [time, power].');
    fclose(fid);
end
end
