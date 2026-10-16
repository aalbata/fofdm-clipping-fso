function fetch_dataset(dataDir)
%FETCH_DATASET Download the power measurements from Zenodo and rebuild the files in dataDir.
%   Source: F. P. Guiomar et al., "Experimental characterization of optical received power
%   in an outdoor free-space optics link," Zenodo, 2020, doi: 10.5281/zenodo.3708496
%   (Creative Commons Attribution 4.0 International). Writes power_sun.csv and power_rain.csv.
url = 'https://zenodo.org/records/3708496/files/MATLAB_FSO-power.zip?download=1';
work = tempname;
mkdir(work);
zipFile = fullfile(work, 'MATLAB_FSO-power.zip');
websave(zipFile, url);
info = dir(zipFile);
if info.bytes ~= 1093380
    error('fetch_dataset:size', 'Unexpected size of the downloaded archive: %d bytes.', info.bytes);
end
unzip(zipFile, work);
write_power_records(fullfile(work, 'MATLAB_FSO-power'), dataDir);
rmdir(work, 's');
end
