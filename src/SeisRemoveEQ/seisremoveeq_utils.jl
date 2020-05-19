"""
    s_whiten!(data::SeisChannel)

Apply spectral whitening to seischannel trace
# Input:
    - `data::SeisChannel`    : SeisChannel from SeisIO
    - `freqmin::Real`: Pass band low corner frequency.
    - `freqmax::Real`: Pass band high corner frequency.
    - `pad::Int`: Number of tapering points outside whitening band.
"""
function s_whiten!(data::SeisChannel, freqmin::Float64, freqmax::Float64;pad::Int=50)

    # compute fft of time series
    FFT = rfft(data.x,1)

    # to use SeisNoise.whiten!() prepare (N, 2) array
    FFTtemp = Array{Complex{Float32}}(undef, length(FFT), 2)
    FFTtemp[:, 1] = FFT

    SeisNoise.whiten!(FFTtemp,freqmin,freqmax,data.fs, data.t[end, 1], pad=pad)

    data.x = irfft(FFTtemp[:,1],data.t[end, 1],1)
    return nothing

end
s_whiten(data::SeisChannel, freqmin::Float64, freqmax::Float64;pad::Int=50) = (U = deepcopy(data);
    s_whiten!(U,freqmin,freqmax,pad=pad);
    return U)
