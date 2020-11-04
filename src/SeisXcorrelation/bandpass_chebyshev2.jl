"""
   bandpass_chebyshev2!(A,freqmin,freqmax,fs,N=5,ripple=40, zerophase=true)
Chebyshev2-Bandpass Filter.
Filter data `A` from `freqmin` to `freqmax` using `N` pole Chebyshev type II filter with `ripple` dB ripple in the stopband.
.
# Arguments
- `A::AbstractArray`: Data to filter
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
- `fs::Real`: Sampling rate in Hz.
- `fs::Int`: Filter corners / order.
- `zerophase::Bool`: If True, apply filter once forwards and once backwards.
This results in twice the filter order but zero phase shift in
the resulting filtered trace.

This function is modified from `bandpass()` in SeisNoise.jl (https://github.com/tclements/SeisNoise.jl/blob/master/src/filter.jl).
"""
function bandpass!(A::AbstractArray{<:AbstractFloat},
                   freqmin::Real, freqmax::Real, fs::Real;
                   N::Int=5, ripple::Int=40, zerophase::Bool=true)
   T = eltype(A)
   fe = T(0.5 * fs)
   low = T(freqmin / fe)
   high = T(freqmax / fe)

   # warn if above Nyquist frequency
   if high - oneunit(high) > -1e-6
       @warn "Selected high corner frequency ($freqmax) of bandpass is at or
       above Nyquist ($fe). Applying a high-pass instead."
       highpass!(A,freqmin,fs,corners=corners,zerophase=zerophase)
       return nothing
   end

   # throw error if low above Nyquist frequency
   if low >= 1
       throw(ArgumentError("Selected low corner frequency is above Nyquist."))
   end

   # create filter
   responsetype = Bandpass(T(freqmin), T(freqmax); fs=fs)
   designmethod = Butterworth(T,corners)

   # use gpu-specific kernel if on the GPU
   if isa(A,AbstractGPUArray)
       gpufilter!(A,responsetype,designmethod)
       return nothing
   end

   # filter if on the CPU
   if zerophase
       A[:,:] .= filtfilt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
   else
       A[:,:] .= filt(digitalfilter(responsetype, designmethod), @view(A[:,:]))
   end

   return nothing
end
bandpass(A::AbstractArray{<:AbstractFloat},freqmin::Real,
         freqmax::Real, fs::Real; corners::Int=4,zerophase::Bool=true) =
         (U = deepcopy(A);bandpass!(U,freqmin,freqmax, fs, corners=corners,
          zerophase=zerophase);return U)
bandpass!(R::RawData,freqmin::Real,freqmax::Real;
          corners::Int=4,zerophase::Bool=true) = (bandpass!(R.x,freqmin,freqmax,
          R.fs,corners=corners,zerophase=zerophase);setfield!(R,:freqmin,Float64(freqmin));
          setfield!(R,:freqmax,Float64(min(freqmax,R.fs/2)));return nothing)
bandpass(R::RawData,freqmin::Real,freqmax::Real;
        corners::Int=4,zerophase::Bool=true) = (U = deepcopy(R);bandpass!(U.x,
        freqmin,freqmax,U.fs,corners=corners,zerophase=zerophase);
        setfield!(U,:freqmin,Float64(freqmin));
        setfield!(U,:freqmax,Float64(min(freqmax,U.fs/2)));return U)
bandpass!(C::CorrData,freqmin::Real,freqmax::Real;
          corners::Int=4,zerophase::Bool=true) = (bandpass!(C.corr,freqmin,freqmax,
          C.fs,corners=corners,zerophase=zerophase);setfield!(C,:freqmin,Float64(freqmin));
          setfield!(C,:freqmax,Float64(min(freqmax,C.fs/2)));return nothing)
bandpass(C::CorrData,freqmin::Real,freqmax::Real;
        corners::Int=4,zerophase::Bool=true) = (U = deepcopy(C);bandpass!(U.corr,
        freqmin,freqmax,U.fs,corners=corners,zerophase=zerophase);
        setfield!(U,:freqmin,Float64(freqmin));
        setfield!(U,:freqmax,Float64(min(freqmax,U.fs/2)));return U)
bandpass!(C::SeisChannel, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=true) = filtfilt!(C,
                   fl=Float64(freqmin),fh=Float64(freqmax),np=corners,rt="Bandpass")
bandpass(C::SeisChannel,freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=true) = (U = deepcopy(C);bandpass!(U,freqmin,freqmax,
                                  corners=corners,zerophase=zerophase);return U)
bandpass!(S::SeisData, freqmin::Real, freqmax::Real;
                   corners::Int=4, zerophase::Bool=true) = filtfilt!(S,
                   fl=Float64(freqmin),fh=Float64(freqmax),np=corners,rt="Bandpass")
bandpass(S::SeisData, freqmin::Real, freqmax::Real; corners::Int=4,
         zerophase::Bool=true) = filtfilt(S,fl=Float64(freqmin),fh=Float64(freqmax),np=corners,rt="Bandpass")
