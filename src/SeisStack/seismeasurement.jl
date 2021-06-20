include("SeisMeasurement/dualstretching.jl")
include("SeisMeasurement/compute_dvvdqq.jl")
include("SeisMeasurement/seisdvv_mwcs.jl")
include("SeisMeasurement/seisdvv_stretching.jl")

"""
    seismeasurement!(C::CorrData, InputDict::OrderedDict)

compute dv/v and dQinv using C.misc["reference"] and C.corr[:,1].

Methods are imported from SeisDvv.jl

# Currently avaliable method:

- `Stretching`      : Canonical time stretching methoc
- `MWCS`            : Moving window cross spectrum
- `WCC`             : Windowed cross-correlation
- `DTW`             : Dynamic time warping
- `DualStretching`  :

NOTE: dQin are obtained only with DualStretching method (dvv is also computed during DualStretching).
"""
function seismeasurement!(C::CorrData, InputDict::OrderedDict)

    !haskey(C.misc, "reference") && (@warn("$(C.name) does not have reference."); return nothing;)

    ref = C.misc["reference"]
    cur = C.corr[:,1]

    measurement_method = InputDict["measurement_method"]

    figdir = InputDict["stretch_debugplot"] ? joinpath(abspath(InputDict["project_outputdir"]), "plots/stack") : ""

	fc = (C.freqmin+C.freqmax)/2.0
	isempty(C.misc["coda_window"]) && return nothing # if coda window is empty, don't perform coda Q measurement

    if lowercase(measurement_method) == "stretching"
		MeasurementDict = seisdvv_stretching(ref, cur, C.misc["timelag"], C.misc["coda_window"],
								C.freqmin,C.freqmax,
								dvmin=-InputDict["dvv_stretching_range"],　dvmax=InputDict["dvv_stretching_range"],
								ntrial_v=InputDict["dvv_stretching_Ntrial"])

    elseif lowercase(measurement_method) == "mwcs"

		MeasurementDict = seisdvv_mwcs(ref,cur,C.freqmin,C.freqmax,C.fs,-C.maxlag,
								InputDict["mwcs_window_length"], InputDict["mwcs_window_step"], InputDict["mwcs_max_dt"],
								InputDict["mwcs_smoothing_half_win"], InputDict["coda_init_factor"], InputDict["max_coda_length"],
								InputDict["min_ballistic_twin"], C.dist*1e3, InputDict["background_vel"])


    elseif lowercase(measurement_method) == "wcc"
        @warn("measurement_method:$(measurement_method) is under developping. skip seismeasurement.")

    elseif lowercase(measurement_method) == "dtw"
        @warn("measurement_method:$(measurement_method) is under developping. skip seismeasurement.")

    elseif lowercase(measurement_method) == "dualstretching"

        fc = (C.freqmin+C.freqmax)/2.0

        MeasurementDict = dualstretching(ref, cur, C.misc["timelag"], fc, C.misc["coda_window"],
                                # parameters for dvv stretching
                                dvmin=-InputDict["dvv_stretching_range"],
                                dvmax=InputDict["dvv_stretching_range"],
                                ntrial_v=InputDict["dvv_stretching_Ntrial"],
                                # parameters for coda-Q stretching
                                dQcinvmin=-InputDict["dQc_stretching_range"],
                                dQcinvmax=InputDict["dQc_stretching_range"],
                                ntrial_q=InputDict["dQc_stretching_Ntrial"],
                                dAAmin=-InputDict["dAA_stretching_range"],
                                dAAmax=+InputDict["dAA_stretching_range"],
                                ntrial_A=InputDict["dAA_stretching_Ntrial"],
                                smoothing_window_len=InputDict["smoothing_window_len"],
                                # parameters for distance and debug figure plot
                                dist_method=InputDict["stretch_distmethod"],
                                figdir=figdir,
                                figname=C.name*"--"*string(C.misc["stack_centraltime"])*"--"*join([C.freqmin, C.freqmax], "-"),
                                fillbox =C.misc["fillbox"])

    elseif lowercase(measurement_method) == "compute_dvvdqq"

		MeasurementDict = compute_dvvdqq(ref, cur, C.misc["timelag"], fc, C.misc["coda_window"],
								C.freqmin,C.freqmax,
								dvmin=-InputDict["dvv_stretching_range"],
								dvmax=InputDict["dvv_stretching_range"],
								ntrial_v=InputDict["dvv_stretching_Ntrial"],
		                        geometrical_spreading_α=InputDict["geometricalspreading_α"],
		                        coda_smooth_window=InputDict["smoothing_window_len"],
		                        figdir=figdir,
		                        figname=C.name*"--"*string(C.misc["stack_centraltime"])*"--"*join([C.freqmin, C.freqmax], "-"),
		                        fillbox=C.misc["fillbox"])

    else
        error("measurement_method:$(measurement_method) is not available.")
    end

    # append measurement results to C.misc
    for key in keys(MeasurementDict)
        C.misc[key] = MeasurementDict[key]
    end

    C.misc["measurement_method"] = measurement_method

    return nothing
end
