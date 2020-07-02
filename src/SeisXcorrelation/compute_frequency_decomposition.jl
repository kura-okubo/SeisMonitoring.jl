"""
    compute_frequency_decomposition(C::CorrData, freqency_band::Array{Float64,1}; cc_bpfilt_method::String="Butterworth",
                                 dj::Float64 = 1 / 12, α0::AbstractFloat = 0.0, αmax::AbstractFloat = 0.25)

compute frequency decomposition of cross-correlation function in corrdata

# Argument
- `C::CoreData` : CorrData contains broadband frequency contents.
- `cc_bpfilt_method::String` : Bandpass filtering method. "Butterworth" or "Wavelet"
- `dj::AbstractFloat`: Spacing between discrete scales. Default value is 1/12.
- `α0::AbstractFloat=0.0`: Lowest tapering fraction for frequency adaptive tapering.
- `αmax::AbstractFloat=0.25`: Highest tapering fraction for frequency adaptive tapering.

# Return
- `C_all::Array{CorrData, 1}`: CorrData contains narrow frequency band
- `freqband::Array{Array{Float64,1},1}` : frequency band used in the C_all.
"""
function compute_frequency_decomposition(
   C_origin::CorrData,
   freqency_band::Array{Float64,1};
   cc_bpfilt_method::String = "ButterWorth",
   dj::AbstractFloat = 1 / 12,
   α0::Float64 = 0.0,
   αmax::Float64 = 0.25
)

   # compute freqband
   Nfreqband = length(freqency_band) - 1
   freqband = map(i -> [freqency_band[i], freqency_band[i+1]], 1:Nfreqband)
   dt = 1.0/C_origin.fs
   C_all = CorrData[]
   C_broadband = deepcopy(C_origin)

   # tapering is applied both before and after cwt filtering to avoid edge effect
   SeisNoise.taper!(C_broadband, max_percentage = 0.1, max_length = 60.0)

   #---ButterWorth---#
   if lowercase(cc_bpfilt_method) == "butterworth"
      for fb in freqband
         freqmin, freqmax = fb
         # apply bandpass filter
         Ctemp = bandpass(C_broadband, freqmin, freqmax, corners = 4)
         # apply frequenc_dependent tapering
         taper_for_freqdecomposition!(Ctemp, freqency_band, α0, αmax)
         Ctemp.misc["freq_decomposition_method"] = cc_bpfilt_method
         push!(C_all, Ctemp)
      end

   #---Wavelet Transform---#
   elseif lowercase(cc_bpfilt_method) == "wavelet"
      # compute cwt and icwt within each frequency band
      T, N = size(C_broadband.corr)

      # 3D corr array to temporally store the all traces
      tr_freq_reconstructed = zeros(Float32, T, N, Nfreqband)

      for traceid = 1:N
         tr = @view C_broadband.corr[:, traceid]
         W, sj, freqs, coi = SeisDvv.cwt(tr,dt,C_broadband.freqmin,C_broadband.freqmax,dj = dj)
         for (ifb, fb) in enumerate(freqband)
            freqmin, freqmax = fb
            # find all index within frequency band
            inds = findall(f -> freqmin <= f <= freqmax, freqs)
            for ind in inds
               inv_W = SeisDvv.icwt(W[:, ind], sj[ind], dt, dj = dj)
               any(isnan.(inv_W)) && (println("nan found in icwt."); println(inv_W))
               #inv_W is reconstructed trace within freqband in time domain
               tr_freq_reconstructed[:, traceid, ifb] = inv_W
            end
         end
      end

      # reshape traces into CorrData
      for (ifb, fb) in enumerate(freqband)
         # copy all metadata
         Ctemp = deepcopy(C_broadband)
         freqmin, freqmax = fb
         Ctemp.freqmin = freqmin
         Ctemp.freqmax = freqmax
         Ctemp.corr = tr_freq_reconstructed[:, :, ifb]
         # apply frequenc_dependent tapering
         taper_for_freqdecomposition!(Ctemp, freqency_band, α0, αmax)
         Ctemp.misc["freq_decomposition_method"] = cc_bpfilt_method
         push!(C_all, Ctemp)
      end

   else
      error("cc_bpfilt_method $(cc_bpfilt_method) is not available. use ButterWorth or Wavelet")
   end

   return (C_all, freqband)
end

"""
   taper_for_freqdecomposition!(C::CorrData, α0::AbstractFloat, αmax::AbstractFloat)

compute frequency adaptive tapering.

```
          /------------\\
         /              \\
  ------/                \\-------
  margin tukeywindow (α=0.1) margin
```
"""
function taper_for_freqdecomposition!(C::CorrData, freqency_band::Array{Float64,1},
                                    α0::AbstractFloat, αmax::AbstractFloat)

   T, N = size(C.corr)
   #adaptive taper percentage due to stronger edge effect at high frequency
   max_percentage = ((αmax-α0)/(freqency_band[end]-freqency_band[2]))*(C.freqmax-freqency_band[2]) + α0
   taperwindow = zeros(T)
   margin = round(Int, 0.5*max_percentage*T)
   tukeywindow = DSP.tukey(T-2*margin, 0.1)
   taperwindow[margin+1:T - margin] = tukeywindow

   for i in 1:N
      C.corr[:, i] .*= taperwindow
   end
end
#
# using  SeisIO, SeisNoise, JLD2, DSP, Dates, DataStructures, SeisDvv
#
# c1 = jldopen("corr_test.jld2", "r")
# C  = c1["test"]
# close(c1)
#
# freqency_band = [0.1, 0.2, 0.5, 0.9]
#
# C_all, freqband = compute_frequency_decomposition(C, freqency_band, cc_bpfilt_method="Wavelet")
#
# using Plots
#
# plot(bg=:white)
#
# for i = 1:length(C_all)
#     C1 = C_all[i]
#     tvec = collect(-C.maxlag:1/C1.fs:C.maxlag)
#     plotid = 10
#     normamp = maximum(abs.(C1.corr[:, plotid]))
#     freqlabel = join([string(freqband[i][1]), string(freqband[i][2])], "-")
#     plot!(tvec, C1.corr[:, plotid]./normamp, label=freqlabel )
# end
# xlabel!("Time lag [s]")
