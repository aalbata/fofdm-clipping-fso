# Clipping and subband filtering for OFDM over a measured free-space optical link

MATLAB code, data, and results for the paper

> S. Alghadi, S. Alkaabi, E. Abd-Elrady, W. Zgallai, and A. Al Bataineh, "Performance Tradeoffs of Clipping
> and Sub-Band Filtering in OFDM over Outdoor Free-Space Optical Channels," accepted at the International
> Near East Conference on Electrical and Energy Systems (NEECON 2026), Dubai, United Arab Emirates,
> 21–23 October 2026.

The code simulates four waveform configurations over received-power records of a 55 m outdoor
free-space optical link: conventional orthogonal frequency-division multiplexing (OFDM), filtered
OFDM (F-OFDM), clipped OFDM, and joint
clipping-assisted F-OFDM, in which a subband filter is followed by an envelope limiter. It reports
the bit error rate (BER), the signal-to-noise ratio (SNR) required for BER targets under two
transmitter models, and the spectrum and peaks of the transmitted stream. The figures of the paper
are produced by the scripts in this repository, and the values of its result tables are those of
`results/clipping_reference.csv` and `results/summary.csv`.

## Requirements

MATLAB R2026b. No toolboxes are required. The results in this repository were produced with
MATLAB R2026b on Windows 11 (Intel Core Ultra 7 155U). The figure export uses options of
`exportgraphics` that are available from R2025a.

## Reproducing the results

Run from the repository folder:

```matlab
reproduce            % tables and figures from the shipped simulation results (less than a minute)
reproduce('full')    % reruns the simulation from the power records in data/ (about 11 minutes)
reproduce('raw')     % downloads the dataset from Zenodo, rebuilds data/, then runs as 'full'
```

The checks of the building blocks are run with

```matlab
cd tests
run_tests
```

The simulation uses fixed random seeds. A rerun on the same machine and MATLAB release reproduces
every file in `results/` byte for byte. The figure files contain creation dates and are therefore
not byte-identical between runs.

## Simulation in brief

| Item | Setting |
|---|---|
| Waveform | 4096-point fast Fourier transform (FFT), 3072 data subcarriers, 15 kHz spacing, 410-sample cyclic prefix, uncoded Gray-mapped 16-ary quadrature amplitude modulation (16-QAM) |
| Subband filter | 512 taps, sinc with Chebyshev window (75 dB sidelobe attenuation), tone offset 2.5 subcarriers |
| Limiter | Envelope limiter at four times the sampling rate; clipping ratio = threshold / root-mean-square (RMS) amplitude of the limiter input, 0.2 to 3.0 in steps of 0.1 |
| Channel | One gain per OFDM block from a power record normalized to unit mean; additive white Gaussian noise |
| Receiver | Matched subband filter (filtered configurations), zero-forcing equalizer with known response, hard decisions |
| Average-power model | Unit average transmit power; SNR = 1 / noise variance |
| Peak-limited model | Unit peak amplitude; peak SNR (PSNR) = 1 / noise variance = SNR x peak-to-average power ratio of the limiter output; OFDM and F-OFDM operate at the clipping ratio 3 |
| Runs | Five runs per record that together use every record sample once; 22.1 million (sunny) and 33.2 million (rainy) bits per BER value |

## Results

Required SNR (average-power model) and required PSNR (peak-limited model) in dB, and guard-band
level of the transmitted stream. The clipped configurations use the clipping ratio 1.7, which
minimizes the required PSNR at a BER of 1e-3.

| Model | Configuration | Clipping ratio | BER 1e-2, sunny | BER 1e-2, rainy | BER 1e-3, sunny | BER 1e-3, rainy | Guard band (dB) |
|---|---|---|---|---|---|---|---|
| Average power | OFDM | none | 12.67 | 12.86 | 15.32 | 15.68 | -36 |
| Average power | F-OFDM | none | 12.67 | 12.86 | 15.33 | 15.69 | -132 |
| Average power | Clipped OFDM | 1.7 | 13.12 | 13.31 | 16.16 | 16.49 | -26 |
| Average power | Joint clipping-assisted F-OFDM | 1.7 | 13.12 | 13.31 | 16.18 | 16.50 | -26 |
| Peak-limited | OFDM | 3.0 | 22.21 | 22.40 | 24.87 | 25.22 | -36 |
| Peak-limited | F-OFDM | 3.0 | 22.21 | 22.41 | 24.88 | 25.24 | -57 |
| Peak-limited | Clipped OFDM | 1.7 | 17.98 | 18.17 | 21.02 | 21.34 | -26 |
| Peak-limited | Joint clipping-assisted F-OFDM | 1.7 | 17.98 | 18.17 | 21.04 | 21.36 | -26 |

Column definitions:

- Model: transmitter model. Average power: all configurations have the same average transmit power
  and the value is the required SNR. Peak-limited: all configurations have the same peak amplitude
  and the value is the required PSNR.
- Clipping ratio: clipping threshold divided by the RMS amplitude of the limiter input; "none" means
  no limiter.
- BER 1e-2 and BER 1e-3, sunny and rainy: SNR or PSNR in dB at which the pooled BER of the sunny or
  rainy record crosses the target, interpolated between the simulated SNR values (1 dB grid).
- Guard band: mean power spectral density over frequencies of 24 MHz and above (up to half the
  sampling rate, 30.72 MHz) relative to the mean in-band density, in dB, for the stream of that row.

Peak of the transmitted stream after interpolation to 16 times the sampling rate, relative to the
clipping threshold, at the clipping ratio 1.7 (`results/transmit_stats.csv`, first six frames,
360 blocks): for clipped OFDM and joint clipping-assisted F-OFDM, a median of
0.34 dB and a maximum of 0.6 dB
with the limiter at four times the sampling rate, against medians of
2.76 dB and 2.75 dB with the
limiter at the sampling rate. The guard-band levels of the streams without limiter and with the
clipping ratio 3 and the maxima of the peak excess are statistical estimates whose uncertainty is of
the order of 0.1 dB to 0.3 dB; the guard-band levels in the table are therefore given in whole
decibels.

## Repository layout

| Path | Content |
|---|---|
| `reproduce.m` | Entry point: analysis, tables, and figures; optionally the simulation and the data download |
| `src/run_simulation.m` | BER simulation of all configurations (`ber.csv`, `waveform_stats.csv`, `channel_stats.csv`) |
| `src/analyze_results.m` | Required SNR, selected clipping ratio, clipping reference, closed-form reference, filter response |
| `src/transmit_spectrum.m` | Power spectral density, guard-band level, power beyond half the sampling rate, and peak excess of the transmitted stream |
| `src/write_summary.m` | Summary table and its differences |
| `src/plot_*.m` | Figures of the paper |
| `src/` (other files) | Waveform, filter, limiter, receiver, and helper functions |
| `tests/run_tests.m` | Checks of the building blocks and of the analysis |
| `data/` | Power records used by the simulation |
| `results/` | Simulation results and derived tables (CSV) |
| `figures/` | Figures as PDF and PNG |

## Result files

All files are comma-separated with a header row. `filter` is 0 without and 1 with the subband
filter. `gamma` is the clipping ratio; `Inf` means no limiter. `NaN` marks a value that does not
exist: a BER target that is not reached, or a quantity that is not defined for that row.

| File | Columns |
|---|---|
| `ber.csv` | `record`, `filter`, `gamma`, `snr_db`, `seed` (run number 1 to 5; the generator seed of a run is 1000 times the record index plus the run number), `bit_errors`, `bits`: bit errors and compared bits of one run at one SNR (average-power model) |
| `waveform_stats.csv` | `record`, `seed`, `filter`, `gamma`, `sigma2` (power at the limiter input), `alpha_ref` and `alpha_meas` (attenuation factor of the limiter, Gaussian reference and measured), `pc_ref` and `p_meas` (limiter output power, reference and measured), `e_x` and `e_err` (energy of the transmitted symbols and of the noise-free symbol error) |
| `channel_stats.csv` | `record`, `n_samples`, `sampling_time_s`, `duration_min` (first to last sample), `power_min_dbm`, `power_max_dbm`, `power_mean_dbm` (mean of the linear power), `gain_min_db`, `gain_max_db` (normalized power gain), `scintillation_index` (variance of the normalized power gain) |
| `required_snr.csv` | `record`, `filter`, `gamma`, `ber_target`, `snr_req_db` (pooled over the runs), `snr_req_seed_min_db`, `snr_req_seed_max_db` (range over the runs), `papr_db` (peak-to-average power ratio of the limiter output), `psnr_req_db` (= `snr_req_db` + `papr_db`) |
| `selected_gamma.csv` | `record`, `filter`, `gamma` (clipping ratio with the lowest required PSNR at a BER of 1e-3), `psnr_req_db` |
| `clipping_reference.csv` | `gamma`, `alpha`, `output_power`, `distortion`, `sdr_db`, `papr_db`, `clip_probability` (Gaussian reference of the limiter, powers relative to the input power), `sdr_inband_ofdm_db`, `sdr_inband_fofdm_db` (simulated signal-to-distortion ratio on the data subcarriers without and with the subband filter) |
| `validation.csv` | `record`, `snr_db`, `bit_errors`, `bits`, `ber_sim` (OFDM without limiter), `ber_theory` (closed form for 16-QAM averaged over the gains of the record) |
| `closed_form_reference.csv` | `ber_target`, `channel` (`awgn`: unit gain; `sun`, `rain`: gains of the record), `snr_req_db` (SNR at which the closed form equals the target) |
| `filter_response.csv` | `k` (subcarrier index), `magnitude_db` (magnitude response of the subband filter) |
| `psd.csv` | `freq_mhz` and the power spectral density in dB relative to the mean in-band density (mean over frequencies up to 20 MHz) for six streams: OFDM and F-OFDM without limiter (`ofdm_db`, `fofdm_db`), with the clipping ratio 3 (`ofdm_gamma0_db`, `fofdm_gamma0_db`), and with the selected clipping ratio (`clipped_ofdm_db`, `joint_db`) |
| `transmit_stats.csv` | `configuration`, `filter`, `gamma`, `limiter_rate` and `stream_rate` (multiples of the sampling rate at the limiter and of the analyzed stream), `guard_band_db`, `outside_band_db` (power beyond half the sampling rate relative to the total power), `peak_excess_median_db`, `peak_excess_max_db` (peak of each block after interpolation to 16 times the sampling rate, relative to the clipping threshold). The streams ending in `_rate1` use a limiter at the sampling rate, and those ending in `_lowpass` restrict the output of the oversampled limiter to the sampling bandwidth |
| `summary.csv` | `model` (`avg`: average power, required SNR; `peak`: peak-limited, required PSNR), `configuration`, `filter`, `gamma`, `req_<record>_<target>_db` (`<target>` is `1e2` or `1e3` for the BER targets 1e-2 and 1e-3), `guard_band_db` (whole decibels) |
| `summary_differences.csv` | `quantity`, `ber_target`, `record`, `difference_db` (difference of two-decimal values: entries of `summary.csv`, an entry and the closed form of the same record, or the closed form of the record and the closed form with unit gain), `seed_min_db`, `seed_max_db` (range of the difference over the runs, from the per-run values before rounding; `NaN` for differences between closed-form values). `joint_minus_clipped_ofdm_peak` equals `joint_minus_clipped_ofdm_avg` apart from the rounding of the entries |

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
