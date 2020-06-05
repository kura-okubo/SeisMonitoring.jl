using SeisNoise, SeisIO, Statistics, DSP, FFTW, Plots
"""
    smooth_withfiltfilt(A::AbstractArray; window_len::Int=7, window::Symbol=:bartlett)

Apply DSP.filtfilt() to smooth the waveform.

# Arguments

- "A::AbstractArray": periodogram of spectrum
- "window_len::Int=7": window length
- "window::Symbol=:bartlett": window type: (see https://juliadsp.github.io/DSP.jl/stable/windows/)
"""
function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end


#NOTE:========================================================================#
# spectrum_coherence and spectrum_deconvolution are defined below to replace
# the smoothing function from SeisNoise.smooth() to SeisMonitoring.smooth_withfiltfilt()
#=============================================================================#

"""
  spectrum_coherence!(F,window_len, water_level)
Apply coherence method to FFTData `F`. Where,
``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_A(ω) ∣  ∣ u_B(ω) ∣}``
# Arguments
- `F::FFTData`: FFTData object of fft'd ambient noise data.
- `window_len::Int`: Number of window length to smooth spectrum.
- `water_level::AbstractFloat`: Regularization parameter for spectral smoothing.
                                0.01 is a common value [Mehta, 2007].
"""
function spectrum_coherence!(F::FFTData, window_len::Int, water_level::AbstractFloat=0.0)

    #smoothF = smooth(abs.(F.fft),half_win)
    smoothF = smooth_withfiltfilt(abs.(F.fft), window_len=window_len, window=:rect)
    reg = water_level .* mean(abs.(F.fft),dims=1)
    smoothF .+= reg
    F.fft ./= smoothF
end
spectrum_coherence(F::FFTData,window_len::Int, water_level::AbstractFloat=0) =
          (U = deepcopy(F); spectrum_coherence!(U,window_len,water_level);
          return U)

"""
  spectrum_deconvolution!(F,window_len, water_level)
Apply deconvolution method to FFTData `F`. Where,
``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_B(ω) ∣^2}``
# Arguments
- `F::FFTData`: FFTData object of fft'd ambient noise data.
- `window_len::Int`: Number of window length to smooth spectrum.
- `water_level::AbstractFloat`: Regularization parameter for spectral smoothing.
                              0.01 is a common value [Mehta, 2007].
"""
function spectrum_deconvolution!(F::FFTData, window_len::Int, water_level::AbstractFloat=0)
    # smoothF = smooth(abs.(F.fft).^2,half_win)
    smoothF = smooth_withfiltfilt(abs.(F.fft).^2, window_len=window_len, window=:rect) #NOTE: squared fft.
    reg = water_level .* mean(abs.(F.fft).^2,dims=1)
    smoothF .+= reg
    F.fft ./= smoothF
end
spectrum_deconvolution(F::FFTData,window_len::Int, water_level::AbstractFloat=0) =
              (U = deepcopy(F);spectrum_deconvolution!(U,window_len,water_level);
              return U)
