function Z = equalized_symbols(p, ctx, fr, xs)
%EQUALIZED_SYMBOLS Noise-free FFT outputs of the counted blocks divided by the known linear response.
%   fr is the frame prepared by PREPARE_FRAME for the filter setting of the transmit chain,
%   and xs the transmitted waveform within the sampling bandwidth returned by
%   TRANSMIT_STREAM. The channel multiplies each sample by the amplitude gain of its block;
%   in the filtered chains the receive filter is matched to the subband filter. Z is
%   Nd-by-M: the FFT outputs divided by the response of the channel gain, the filter
%   cascade, and the timing advance. Dividing Z by the attenuation factor of the chain gives
%   the equalized symbols. The transmit scale factor cancels and is therefore omitted.
y = fr.an .* xs;
if fr.filt
    y = filter_stream(y, ctx.g);
end
Z = receive_frame(p, y, fr.Dc) ./ fr.eqz;
end
