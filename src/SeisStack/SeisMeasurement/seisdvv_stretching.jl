using Interpolations, Statistics, Distances, StatsBase, ColorSchemes, Plots
using SeisMonitoring: smooth_withfiltfilt
"""
    seisdvv_stretching(ref::AbstractArray,cur::AbstractArray, time::AbstractArray, window::AbstractArray,
                       fmin::Float64,fmax::Float64,
                       dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500)

    Stretching method using SeisDvv.stretching.
"""
function seisdvv_stretching(ref::AbstractArray,cur::AbstractArray, time::AbstractArray, window::AbstractArray,
                   fmin::Float64,fmax::Float64;
                   dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500)


    MeasurementDict = Dict()

    # This scripts are taken from Nextjournal written by Jared Bryan
    # https://nextjournal.com/jtbryan/comparing-approaches-to-measuring-seismic-
    # phase-variations-in-the-time-frequency-and-wavelet-domains?token=SH2k5KRceh9VFTZmLFG5yJ

    # NOTE: Since if the correlation coefficient between ref and cur trace is not enough,
    # it causes error due to bad fitting using glm.
    if cor(ref, cur) < 0.0
        # skip this pair without appending anything
        return MeasurementDict
    end

    dvv_ts, cc_ts, cdp_Ts, eps_ts, err_ts, allC_ts = SeisDvv.stretching(ref, cur, time,
                        window, fmin, fmax, dvmin=dvmin, dvmax=dvmax, ntrial=ntrial_v);

    MeasurementDict["dvv_ts"]     = dvv_ts
    MeasurementDict["cc_ts"]      = cc_ts
    MeasurementDict["cdp_Ts"]   = cdp_Ts
    MeasurementDict["eps_ts"]   = eps_ts
    MeasurementDict["err_ts"]     = err_ts
    MeasurementDict["allC_ts"] = allC_ts

    return MeasurementDict
end
