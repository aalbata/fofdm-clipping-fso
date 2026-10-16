# Clipping and subband filtering for OFDM over a measured free-space optical channel

MATLAB code, data, and results for the paper

> S. Alghadi, S. Alkaabi, E. Abd-Elrady, W. Zgallai, and A. Al Bataineh, "Power and spectral tradeoffs of
> clipping and subband filtering for OFDM in measurement-informed free-space optical channels."

The code compares orthogonal frequency-division multiplexing (OFDM) and filtered OFDM (F-OFDM)
without a limiter with four transmit chains that add an envelope limiter: clipped OFDM, a subband
filter followed by the limiter (filter-then-clip), the limiter followed by the subband filter
(clip-then-filter), and repeated clipping and frequency-domain filtering of each OFDM symbol (RCF).
The channel gains come from measured received-power records of a 55 m outdoor free-space optical
link in sunny and rainy weather. The code reports the bit error rate (BER), the signal-to-noise ratio
(SNR) and peak SNR (PSNR) required for BER targets, the peak-to-average power ratio (PAPR) of the
transmitted waveform, and its power spectral density. Figs. 2 to 6 of the paper are drawn by the
scripts in this repository from the files in `results/`, and Tables 2 to 4 report values of these
files. Fig. 1 is drawn by `src/plot_block_diagram.m`, and Table 1 lists the settings of
`src/sim_params.m` and the record lengths in `results/channel_stats.csv`.

## Requirements

MATLAB R2026b. No toolboxes are required. The results in this repository were produced with
MATLAB R2026b on Windows 11 (Intel Core Ultra 7 155U). The figure export uses options of
`exportgraphics` that are available from R2025a.

## Reproducing the results

Run from the repository folder:

```matlab
reproduce            % tables and figures from the shipped results (about 2 minutes)
reproduce('full')    % characterization and simulation from the power records in data/ (about 55 minutes)
reproduce('raw')     % downloads the dataset from Zenodo, rebuilds data/, then runs as 'full'
```

The checks of the building blocks are run with

```matlab
cd tests
run_tests
```

The characterization and the simulation use fixed random seeds. A rerun on the same machine and
MATLAB release reproduces every file in `results/` byte for byte. The figure files contain creation
dates and are therefore not byte-identical between runs.

## Simulation in brief

| Item | Setting |
|---|---|
| Waveform | 4096-point fast Fourier transform, 3072 data subcarriers, 15 kHz spacing, 410-sample cyclic prefix, uncoded Gray-mapped 16-ary quadrature amplitude modulation |
| Subband filter | 512 taps, sinc with Chebyshev window (75 dB sidelobe attenuation), half amplitude at subcarrier index ±1539, 3 subcarrier spacings beyond the outermost data subcarrier |
| Limiter | Envelope limiter at four times the sampling rate; clipping ratio = threshold / root-mean-square amplitude of the limiter input, 0.2 to 3.0 in steps of 0.1 |
| RCF | 4 iterations of clipping at four times oversampling and removal of all components outside the data subcarriers, per OFDM symbol |
| Peak power | Instantaneous power exceeded with probability 1e-4 by the transmitted waveform reconstructed at 16 times the sampling rate |
| Calibration | Average transmit power, attenuation factor, peak power, guard-band level, and power beyond half the sampling rate of every chain and clipping ratio from 12 frames with a seed that no simulation run uses |
| Spectrum | Power spectral density of Fig. 5 from 60 frames; Welch estimate with segments of 8192 sampling periods (32768 samples for clipped OFDM and filter-then-clip, which are transmitted at four times the sampling rate), Kaiser window (shape parameter 30), hop of a quarter segment |
| Channel | One amplitude gain per OFDM block from a power record normalized to unit mean; additive white Gaussian noise; coherent detection |
| Receiver | Matched subband filter for the filtered chains, zero-forcing equalizer with known response, hard decisions |
| Average-power model | Unit average transmit power; SNR = 1 / noise variance |
| Peak-limited model | Unit peak power; PSNR = SNR x PAPR |
| Runs | Five runs per record that together use every record sample once; 22.1 million (sunny) and 33.2 million (rainy) bits per BER value; SNR from 0 dB to 36 dB in steps of 0.25 dB |

## Results

Required SNR (average-power model) and required PSNR (peak-limited model) in dB for BER targets of
1e-2 and 1e-3, PAPR of the transmitted waveform in dB, and guard-band level in dB
(`results/summary.csv`). Each chain uses the clipping ratio with the lowest required PSNR at a BER of
1e-3 for the sunny record (`results/selected_gamma.csv`).

| Model | Configuration | Clipping ratio | BER 1e-2, sunny | BER 1e-2, rainy | BER 1e-3, sunny | BER 1e-3, rainy | PAPR | Guard-band level |
|---|---|---|---|---|---|---|---|---|
| Average power | OFDM | none | 12.67 | 12.86 | 15.32 | 15.67 | 9.65 | -36 |
| Average power | F-OFDM | none | 12.67 | 12.86 | 15.33 | 15.69 | 9.66 | -132 |
| Average power | Clipped OFDM | 1.7 | 13.12 | 13.31 | 16.16 | 16.48 | 5.14 | -26 |
| Average power | Filter-then-clip | 1.7 | 13.12 | 13.31 | 16.18 | 16.50 | 5.14 | -26 |
| Average power | Clip-then-filter | 1.8 | 12.95 | 13.14 | 15.85 | 16.19 | 7.18 | -125 |
| Average power | Repeated clipping and filtering (RCF) | 1.9 | 13.14 | 13.33 | 16.22 | 16.54 | 6.32 | -36 |
| Peak-limited | OFDM | none | 22.32 | 22.51 | 24.97 | 25.33 | 9.65 | -36 |
| Peak-limited | F-OFDM | none | 22.33 | 22.52 | 24.99 | 25.35 | 9.66 | -132 |
| Peak-limited | Clipped OFDM | 1.7 | 18.25 | 18.44 | 21.30 | 21.62 | 5.14 | -26 |
| Peak-limited | Filter-then-clip | 1.7 | 18.25 | 18.45 | 21.31 | 21.64 | 5.14 | -26 |
| Peak-limited | Clip-then-filter | 1.8 | 20.12 | 20.32 | 23.03 | 23.37 | 7.18 | -125 |
| Peak-limited | Repeated clipping and filtering (RCF) | 1.9 | 19.46 | 19.66 | 22.54 | 22.87 | 6.32 | -36 |

Column definitions:

- Model: average power, all configurations have the same average transmit power and the value is
  the required SNR; peak-limited, all configurations have the same peak power and the value is the
  required PSNR.
- Clipping ratio: clipping threshold divided by the root-mean-square amplitude of the limiter input;
  "none" means no limiter.
- BER 1e-2 and BER 1e-3, sunny and rainy: SNR or PSNR in dB at which the BER pooled over the five runs
  crosses the target, interpolated on a parabola in log BER through three adjacent grid points.
- PAPR: instantaneous power exceeded with probability 1e-4 by the transmitted waveform reconstructed
  at 16 times the sampling rate, relative to its average power, in dB.
- Guard-band level: mean power spectral density over 24 MHz <= |f| <= 30.72 MHz (both band edges) relative
  to the mean density over |f| <= 20 MHz, in dB, from the 12 calibration frames.

## Where the paper's items come from

| Item of the paper | Script | Results file |
|---|---|---|
| Table 1 | `src/sim_params.m`, `src/run_simulation.m` | `channel_stats.csv` |
| Fig. 1 | `src/plot_block_diagram.m` | none |
| Fig. 2 | `src/plot_clipping_ratio_sweep.m` | `required_snr.csv`, `selected_gamma.csv` |
| Fig. 3 | `src/plot_peak_ccdf.m` | `ccdf.csv`, `selected_gamma.csv` |
| Fig. 4 | `src/plot_ber_vs_snr.m` | `ber.csv`, `calibration.csv`, `validation.csv`, `selected_gamma.csv` |
| Fig. 5 | `src/plot_transmit_psd.m` | `psd.csv` |
| Fig. 6 | `src/plot_tradeoff.m` | `required_snr.csv`, `calibration.csv`, `selected_gamma.csv` |
| Table 2 | `src/analyze_results.m` | `clipping_reference.csv` |
| Table 3 | `src/write_summary.m` | `summary.csv`, `summary_differences.csv` |
| Table 4 | `src/write_summary.m` | `peak_gain_sensitivity.csv` |
| Filter response, channel statistics, and calibration (Sections 2.1 and 2.5) | `src/analyze_results.m`, `src/run_simulation.m`, `src/characterize_transmitters.m` | `filter_response.csv`, `channel_stats.csv`, `calibration.csv` |
| Validation and checks (Sections 3.1 and 3.2) | `src/analyze_results.m`, `src/implementation_checks.m`, `src/write_summary.m` | `validation.csv`, `closed_form_reference.csv`, `implementation_checks.csv`, `summary_differences.csv` |
| Sensitivity (Section 3.5) | `src/characterize_transmitters.m`, `src/transmit_spectrum.m`, `src/analyze_results.m`, `src/write_summary.m` | `peak_convergence.csv`, `calibration.csv`, `psd_check.csv`, `required_snr.csv`, `peak_gain_sensitivity.csv` |

## Repository layout

| Path | Content |
|---|---|
| `reproduce.m` | Entry point: characterization, simulation, analysis, tables, and figures, and optionally the data download |
| `src/sim_params.m` | Constants of the waveform, the transmit chains, and the receiver, the simulation grids, the run indices, and the seeds of the characterization and the bootstrap |
| `src/characterize_transmitters.m` | Calibration of every transmit chain (`calibration.csv`, `peak_convergence.csv`, `ccdf.csv`) |
| `src/run_simulation.m` | BER simulation (`ber.csv`, `channel_stats.csv`) |
| `src/analyze_results.m` | Required SNR and PSNR, selected clipping ratios, reference tables |
| `src/transmit_spectrum.m` | Spectral density of the transmitted waveforms and checks of the estimator |
| `src/implementation_checks.m` | Checks of simplifications of the transmitter and receiver models, of the noise scaling of the receiver, and of the PAPR without the subband filter (`implementation_checks.csv`) |
| `src/write_summary.m` | Summary table, its differences, and the peak-power gains at three exceedance probabilities (`summary.csv`, `summary_differences.csv`, `peak_gain_sensitivity.csv`) |
| `src/plot_*.m` | Figures of the paper |
| `src/` (other files) | Waveform, filter, limiter, transmit chains, receiver, and helper functions |
| `tests/run_tests.m` | Checks of the building blocks and of the analysis |
| `data/` | Power records used by the simulation |
| `results/` | Results and derived tables (CSV) |
| `figures/` | Figures as PDF and PNG |

## Result files

All files are comma-separated with a header row. `chain` is one of `clip`, `filter_clip`,
`clip_filter`, and `rcf`; without a limiter (`gamma` = `Inf`), `clip` and `rcf` are conventional OFDM
and `filter_clip` and `clip_filter` are F-OFDM. `gamma` is the clipping ratio. `NaN` marks a value that
does not exist, such as a BER target that is not reached.

| File | Columns |
|---|---|
| `ber.csv` | `record`, `chain`, `gamma`, `snr_db`, `bits_per_run`, `errors_run1` to `errors_run5`: bit errors of each run at one SNR, average-power model (the generator seed of run i of record r is 1000 r + i, with r = 1 for the sunny and r = 2 for the rainy record) |
| `channel_stats.csv` | `record`, `n_samples`, `sampling_time_s`, `duration_min`, `power_min_dbm`, `power_max_dbm`, `power_mean_dbm` (mean of the linear power), `gain_min_db`, `gain_max_db` (normalized power gain), `scintillation_index` (variance of the normalized power gain) |
| `calibration.csv` | `chain`, `gamma`, `limiter_sigma2` (power at the limiter input), `tx_power` (average transmitted power), `alpha` (attenuation factor), `alpha_gauss` and `power_gauss` (Gaussian reference where it applies), `inband_sdr_db` (signal-to-distortion ratio on the data subcarriers), `papr_1e3_db`, `papr_1e4_db`, `papr_1e5_db` (PAPR at exceedance probabilities 1e-3, 1e-4, 1e-5), `papr_max_db` (largest sample), `papr_1e3_lo_db` to `papr_1e5_hi_db` (95% bootstrap intervals over the calibration frames), `guard_band_db` (guard-band level over the calibration frames), `outside_band_db` (power beyond half the sampling rate relative to the total power) |
| `peak_convergence.csv` | `chain`, `gamma`, `rate` (16 or 32 times the sampling rate), `papr_1e3_db`, `papr_1e4_db`, `papr_1e5_db` from the first four characterization frames |
| `ccdf.csv` | `chain`, `gamma`, `level_db` (instantaneous power relative to the average power), `exceedance` (fraction of samples above that level) |
| `required_snr.csv` | `record`, `chain`, `gamma`, `ber_target`, `snr_req_db` (pooled over the runs), `snr_req_run_min_db`, `snr_req_run_max_db` (range over the runs), `snr_req_1db_grid_db` (from the whole-decibel grid), `snr_req_linear_db` (straight-line interpolation), `papr_db`, `psnr_req_db` (= `snr_req_db` + `papr_db`) |
| `selected_gamma.csv` | `chain`, `gamma` (lowest required PSNR at 1e-3, sunny record), `psnr_req_db`, `near_optimal_min`, `near_optimal_max`, `near_optimal_count` (clipping ratios within 0.1 dB of the minimum), `gamma_rain`, `gamma_sun_other_target`, `gamma_rain_other_target` |
| `clipping_reference.csv` | `gamma`, `alpha`, `output_power`, `distortion`, `sdr_db`, `papr_db`, `clip_probability` (Gaussian reference of the limiter samples, powers relative to the input power), `alpha_clip`, `inband_sdr_clip_db`, `papr_1e4_clip_db` (measured for clipped OFDM) |
| `validation.csv` | `record`, `snr_db`, `bit_errors`, `bits`, `ber_sim` (OFDM without a limiter), `ber_theory` (closed form for 16-QAM averaged over the gains of the record) |
| `closed_form_reference.csv` | `ber_target`, `channel` (`awgn`: unit gain; `sun`, `rain`: gains of the record), `snr_req_db` |
| `filter_response.csv` | `k` (subcarrier index), `magnitude_db` (magnitude response of the subband filter) |
| `psd.csv` | `freq_mhz` and the power spectral density in dB relative to the mean density for \|f\| <= 20 MHz of OFDM and F-OFDM without a limiter (`ofdm_db`, `fofdm_db`) and of each chain at its selected clipping ratio (`clip_db`, `filter_clip_db`, `clip_filter_db`, `rcf_db`), from 60 frames |
| `psd_check.csv` | `stream`, `estimator` (`welch` or `closed_form`), `segment`, `kaiser_beta`, `frames`, `guard_band_db` |
| `implementation_checks.csv` | `check`, `configuration`, `gamma`, `value`: `gain_order` (energy of the noise-free FFT outputs over the energy of their change when the block gains are applied before the restriction to the sampling bandwidth, in dB), `noise_power` (normalized noise power at the FFT output, expected value one), `hard_limit` (energy of the noise-free FFT outputs over the energy of their change when the reconstructed waveform is limited at its peak power, in dB), `clipped_fraction` (fraction of reconstructed samples limited in that check), `lowpass_papr` (PAPR of the clip chain restricted to the sampling bandwidth without a subband filter, in dB) |
| `summary.csv` | `model` (`avg` or `peak`), `configuration`, `chain`, `gamma`, `req_<record>_<target>_db` (`<target>` is `1e2` or `1e3`), `papr_db`, `guard_band_db` (whole decibels), `outside_band_db` |
| `peak_gain_sensitivity.csv` | `chain`, `gamma`, `reference` (unclipped reference chain), `ber_target`, `record`, `gain_papr_1e3_db`, `gain_papr_1e4_db`, `gain_papr_1e5_db` (reduction of the required PSNR relative to the reference, from two-decimal PSNR values, with the peak power at each exceedance probability) |
| `summary_differences.csv` | `quantity`, `ber_target`, `record`, `difference_db` (difference of the two-decimal entries of `summary.csv`, or of an entry and the closed form), `run_min_db`, `run_max_db` (range of the same difference over the individual runs, formed from unrounded values and then rounded, so that `difference_db`, a difference of rounded entries, can lie slightly outside it) |

## Data

The files `data/power_sun.csv` and `data/power_rain.csv` are derived from

> F. P. Guiomar, A. Lorences Riesgo, D. Ranzal, F. Rocco, A. Carena, and P. P. Monteiro,
> "Experimental characterization of optical received power in an outdoor free-space optics link,"
> Zenodo, 2020, doi: 10.5281/zenodo.3708496,

which is distributed under the Creative Commons Attribution 4.0 International license
(https://creativecommons.org/licenses/by/4.0/). Changes made: the variables `power` (received
optical power in dBm) and `samplingTime` (4 s) of the two measurement files were written as text
with the columns `time_s` and `power_dbm`, with the power at the 0.01 dBm resolution of the
measurement. `reproduce('raw')` repeats this conversion from the original archive.

## License

The code is released under the MIT License (see `LICENSE`). The files in `data/` remain under the
Creative Commons Attribution 4.0 International license of the dataset (see `data/README.md`).
