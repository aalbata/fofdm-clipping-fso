function p = sim_params()
%SIM_PARAMS Constants of the waveform, the transmit chains, the receiver, and the grids.
p.N = 4096;                          % FFT size
p.k = [-1536:-1, 1:1536].';          % data subcarrier indices
p.Nd = numel(p.k);                   % number of data subcarriers
p.bins = mod(p.k, p.N) + 1;          % FFT bins of the data subcarriers (1-based)
p.Ncp = 410;                         % cyclic-prefix length (samples)
p.Nb = p.N + p.Ncp;                  % samples per block
p.df = 15e3;                         % subcarrier spacing (Hz)
p.Fs = p.N * p.df;                   % sampling rate (Hz)
p.L = 512;                           % subband filter length (taps)
p.toneOffset = 2.5;                  % excess bandwidth on each side (subcarrier spacings)
p.sidelobeDb = 75;                   % sidelobe attenuation of the Chebyshev window (dB)
p.delta = p.Ncp / 2;                 % timing advance of the FFT window (samples)
p.M = 60;                            % counted blocks per frame
p.B = p.M + 2;                       % blocks per frame, with one guard block at each end
p.os = 4;                            % oversampling factor of the limiter
p.nPad = 2^19;                       % frame length after zero padding (samples at the sampling rate)

% Transmit chains. Without a limiter (clipping ratio Inf), clip and rcf give conventional
% OFDM and filter_clip and clip_filter give filtered OFDM (F-OFDM).
%   clip         limiter on the OFDM stream
%   filter_clip  subband filter, then limiter
%   clip_filter  limiter, restriction to the sampling bandwidth, then subband filter
%   rcf          repeated clipping and frequency-domain filtering of each OFDM symbol
p.chains = {'clip', 'filter_clip', 'clip_filter', 'rcf'};
p.chainFiltered = [false, true, true, false];     % subband filter in the chain and matched receive filter
p.rcfIterations = 4;                 % clipping and filtering iterations of the rcf chain

p.gammaGrid = [(2:30) / 10, Inf];    % clipping ratios; Inf means no limiter
p.snrGrid = 0:0.25:36;               % SNR grid (dB)
p.seeds = 1:5;                       % runs
p.records = {'sun', 'rain'};         % power records
p.berTargets = [1e-2, 1e-3];         % BER targets for the required SNR
p.selectionTarget = 1e-3;            % BER target used to select the clipping ratio
p.nearOptimalDb = 0.1;               % margin of the near-optimal clipping ratios (dB)

% Characterization of the transmitted waveform.
p.peakRate = 16;                     % sampling-rate multiple of the reconstructed waveform
p.peakRateCheck = 32;                % sampling-rate multiple of the convergence check
p.peakProb = 1e-4;                   % exceedance probability that defines the peak power
p.peakProbs = [1e-3, 1e-4, 1e-5];    % exceedance probabilities reported
p.calSeed = 9001;                    % seed of the characterization frames (not used by any run)
p.levelBinDb = 0.01;                 % bin width of the power histograms (dB)
p.checkGamma = [1.0 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.5 3.0 Inf];  % clipping ratios of the rate check
p.checkFrames = 4;                   % characterization frames of the rate check
p.bootSeed = 7;                      % seed of the bootstrap over frames
p.bootDraws = 1000;                  % bootstrap resamples
p.calFrames = 12;                    % frames for the attenuation factor, power, peaks, and guard band
p.psdFrames = 60;                    % frames for the spectral densities of the figures
p.psdSegment = 8192;                 % Welch segment length (sampling periods)
p.kaiserBeta = 30;                   % shape parameter of the Kaiser window
p.inBandHz = 20e6;                   % in-band normalization: |f| <= 20 MHz
p.guardHz = 24e6;                    % guard-band level range: 24 MHz <= |f| <= Fs/2
end
