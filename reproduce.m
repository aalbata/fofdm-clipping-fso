function reproduce(mode)
%REPRODUCE Regenerate the tables and figures of the paper.
%   reproduce           tables and figures from the shipped simulation results
%   reproduce('full')   rerun the simulation from the power records in data/ first
%   reproduce('raw')    download the dataset from Zenodo and rebuild data/ first, then as 'full'
if nargin < 1
    mode = 'fast';
end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
dataDir = fullfile(root, 'data');
resultsDir = fullfile(root, 'results');
figuresDir = fullfile(root, 'figures');

switch mode
    case 'fast'
        simulate = false;
    case 'full'
        simulate = true;
    case 'raw'
        fetch_dataset(dataDir);
        simulate = true;
    otherwise
        error('reproduce:mode', 'Unknown mode "%s". Use ''fast'', ''full'', or ''raw''.', mode);
end

if simulate
    run_simulation(dataDir, resultsDir);
end
analyze_results(dataDir, resultsDir);
if simulate
    transmit_spectrum(resultsDir);
end
write_summary(resultsDir);

plot_block_diagram(figuresDir);
plot_clipping_ratio_sweep(resultsDir, figuresDir);
plot_ber_vs_snr(resultsDir, figuresDir);
plot_transmit_psd(resultsDir, figuresDir);
fprintf('Tables are in %s\nFigures are in %s\n', resultsDir, figuresDir);
end
