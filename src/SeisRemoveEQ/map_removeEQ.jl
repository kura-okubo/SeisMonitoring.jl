include("seisremoveeq_utils.jl")
include("seisremoveeq_kurtosis.jl")
include("seisremoveeq_stalta.jl")
include("seisremoveeq_remove_eqfilt.jl")

"""
    map_removeEQ(station, InputDict::OrderedDict)
    remove earthquake and save it into jld2 file.

"""
function map_removeEQ(stationpath::String, InputDict::OrderedDict)

	fidir, fikey = splitdir(stationpath)

	# process only .seisio file
	if split(fikey, ".")[end] != "seisio"
		return (0, 0)
	end

    println("start process on $(fikey)")

    # open file io at each process
    # fi = jldopen(InputDict["RawData_path"], "r")

    bt_kurtosis = 0.0
    bt_stalta = 0.0

    # for fikey in keys(fi["Waveforms/$(station)"])

	bt_1, bt_2, btsta_1, bt_3 = zeros(4)
    #read SeisChannel data
    # S1_raw = fi["Waveforms/$(station)/$(fikey)"]
	S1_raw = rseis(stationpath)[1]
    S1 = deepcopy(S1_raw)
    # skip if it does not exist
    if isempty(S1.x) || all(isnan.(S1.x))
        return (bt_kurtosis, bt_stalta)
    end

    #NOTE:Transient error on combination between download margin and resampling#
    # S.t may contain index 0 in t[:, 1]. In that case, this causes the error on later part, so discard that data.
    if any(S1.t[:, 1] .== 0)
        @warn("zero index found. discard it.")
        return (bt_kurtosis, bt_stalta)
    end

    #initiate removal flags
    S1.misc["noisesignal"] = fill(true, length(S1.x))

    if InputDict["IsKurtosisRemoval"]
        # compute kurtosis and detect
        bt_1 = @elapsed get_kurtosis!(S1, InputDict)
        bt_2 = @elapsed detect_eq_kurtosis!(S1, InputDict)
    end

    if InputDict["IsSTALTARemoval"]
        # detect earthquake and tremors by STA/LTA
        btsta_1 = @elapsed detect_eq_stalta!(S1, InputDict)
    end

    # apply spectral whitening before remove filtering
    InputDict["IsWhitening"] && s_whiten!(S1,InputDict["freqmin_whiten"], InputDict["freqmax_whiten"])

    bt_3 = @elapsed remove_eqfilt!(S1, InputDict)

    bt_kurtosis += bt_1 + bt_2
    bt_stalta += btsta_1

    #compute removal fraction on this channel
    ns_list = S1.misc["noisesignal"]
    numofremoval = sum(x -> x == false, ns_list, dims = 1)
    S1.misc["removal_fraction"] = float(numofremoval[1]) /
                                  float(length(ns_list))
	# append noise data fraction
	S1.misc["data_fraction"] = get_noisedatafraction(S1.x, zerosignal_minpts=100, eps_α=1e-6)

    # Append raw trace, kurtosis and stalta time traces to S1
    if InputDict["Append_alltraces"]
        S1.misc["raw_trace"] = S1_raw.x
        S1.misc["kurtsis_trace"] = S1.misc["kurtosis"]
        S1.misc["stalta_trace"] = S1.misc["stalta"]
        delete!(S1.misc, "kurtosis")
        delete!(S1.misc, "stalta")
		# NOTE: misc cannot saved with seisio native binary, so save with JLD2.
		fopath = joinpath(InputDict["tmpdir"], fikey * ".jld2")
		fo = jldopen(fopath, true, true, true, IOStream)
		fo["S"] = S1
		JLD2.close(fo)
    else
        # don't keep those traces to reduce data size
        delete!(S1.misc, "kurtosis")
        delete!(S1.misc, "stalta")
        delete!(S1.misc, "noisesignal")

		# write with native format for next process
		temppath = joinpath(InputDict["tmpdir"], fikey)
		wseis(temppath, S1)
    end

	#DEBUG: to avoid too much memory allocation, save each seischannel into seisio file

	# # temppath = joinpath(InputDict["tmpdir"], fikey * ".jld2")
	# temppath = joinpath(InputDict["tmpdir"], fikey * ".seisio")
	#
	# wseis(temppath, S1)

	# jldopen(temppath, "w") do fo
		# fo["S"] = S1
	# end
    # end

    # JLD2.close(fi)

    return (bt_kurtosis, bt_stalta)
end
