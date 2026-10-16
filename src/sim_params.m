function p = sim_params()
%SIM_PARAMS Constants of the waveform, the receiver, and the simulation grids.
p.N = 4096;                          % FFT size
p.k = [-1536:-1, 1:1536].';          % data subcarrier indices
p.Nd = numel(p.k);                   % number of data subcarriers
p.bins = mod(p.k, p.N) + 1;          % FFT bins of the data subcarriers (1-based)
p.Ncp = 410;                         % cyclic-prefix length (samples)
p.Nb = p.N + p.Ncp;                  % samples per block
p.df = 15e3;                         % subcarrier spacing (Hz)
p.Fs = p.N * p.df;                   % sampling rate (Hz)
p.L = 512;                           % subband filter length (taps)
p.toneOffset = 2.5;                  % excess bandwidth on each side (subcarriers)
p.sidelobeDb = 75;                   % sidelobe attenuation of the Chebyshev window (dB)
p.delta = p.Ncp / 2;                 % timing advance of the FFT window (samples)
p.M = 60;                            % counted blocks per frame
p.B = p.M + 2;                       % blocks per frame, with one guard block at each end
p.os = 4;                            % oversampling factor of the limiter
p.nPad = 2^19;                       % frame length after zero padding (samples at the sampling rate)
p.gammaGrid = [(2:30) / 10, Inf];    % clipping ratios; Inf means no limiter
p.gamma0 = 3;                        % clipping ratio without deliberate clipping (peak-limited model)
p.snrGrid = 0:36;                    % SNR grid (dB)
p.seeds = 1:5;                       % independent runs
p.records = {'sun', 'rain'};         % power records
p.berTargets = [1e-2, 1e-3];         % BER targets for the required SNR
p.selectionTarget = 1e-3;            % BER target used to select the clipping ratio
end
