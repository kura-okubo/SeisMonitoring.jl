using Interpolations, Statistics, Distances, StatsBase, ColorSchemes, Plots
using SeisMonitoring: smooth_withfiltfilt
"""
    seisdvv_mwcs(ref::AbstractArray, cur::AbstractArray, t::AbstractArray, fc::Float64,
                        window::AbstractArray;
                        dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500,
                        dQcinvmin::Float64=-0.01, dQcinvmax::Float64=+0.01, ntrial_q::Int=101,
                        dAAmin::Float64=-0.02, dAAmax::Float64=+0.02, ntrial_A::Int=101,
                        dist_method::String="euclidean")

    MWCS method using SeisDvv.mwcs amd SeisDvv.mwcs_dvv.
    Use dynamic time lag.
"""
function seisdvv_mwcs(ref::AbstractArray,cur::AbstractArray,fmin::Float64,
                   fmax::Float64,fs::Float64,tmin::Float64,
                   window_length::Float64,window_step::Float64, max_dt::Float64,
                   smoothing_half_win::Int, coda_init_factor::Real, max_coda_length::Real, min_ballistic_twin::Real,
                   dist::Float64, dtt_v::Float64)

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

    t_axis_mwcs, dt_mwcs, error_mwcs, mcoh_mwcs = SeisDvv.mwcs(ref, cur, fmin,
            fmax, fs, tmin, window_length, window_step, smoothing_half_win);

    # NOTE: As experiment, dtt_width is defined as a factor of "mwcs_window_length"
    # NOTE: As experiment, fix max_dt at 1.0 [s].
    # dtt_width = 6 * window_length
    # dvv_mwcs, dvv_err_mwcs, int_mwcs, int_err_mwcs, dvv0_mwcs, dvv0_err_mwcs = SeisDvv.mwcs_dvv(t_axis_mwcs,
    #  dt_mwcs, error_mwcs, mcoh_mwcs, "dynamic", dist, dtt_v, 0.0, dtt_width, "both", max_dt=1.0);

    dtt_minlag = max(min_ballistic_twin, coda_init_factor * dist / dtt_v) #[s]
    dtt_width = max_coda_length - dtt_minlag
    # dtt_width has been sure to have enough lenght as slice_codawindow() already rejects it if it's too small.
    dvv_mwcs, dvv_err_mwcs, int_mwcs, int_err_mwcs, dvv0_mwcs, dvv0_err_mwcs = SeisDvv.mwcs_dvv(t_axis_mwcs,
        dt_mwcs, error_mwcs, mcoh_mwcs, "static", dist, dtt_v, dtt_minlag, dtt_width, "both", max_dt=max_dt);

     #DEBUG: max_dt is uncertain as findall(x -> abs.(x) .>= max_dt,time_axis) l 244 can be findall(x -> abs.(x) .>= max_dt, dt)

    MeasurementDict["t_axis_mwcs"]  = t_axis_mwcs
    MeasurementDict["dt_mwcs"]      = dt_mwcs
    MeasurementDict["error_mwcs"]   = error_mwcs
    MeasurementDict["mcoh_mwcs"]    = mcoh_mwcs
    MeasurementDict["dvv_mwcs"]     = dvv_mwcs
    MeasurementDict["dvv_err_mwcs"] = dvv_err_mwcs
    MeasurementDict["dvv0_mwcs"]    = dvv0_mwcs
    MeasurementDict["dvv0_err_mwcs"]= dvv0_err_mwcs

    # MeasurementDict["s_alltrace"] = s_alltrace # please avoid to save 3D array for the moment if the potential issue occurs in cluster.
    return MeasurementDict
end
