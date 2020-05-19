include("seisremoveeq_utils.jl")
include("seisremoveeq_kurtosis.jl")
include("seisremoveeq_stalta.jl")
include("seisremoveeq_remove_eqfilt.jl")

"""
    map_removeEQ(station, InputDict::OrderedDict)
    remove earthquake and save it into jld2 file.

"""
function map_removeEQ(station::String, InputDict::OrderedDict)

    println("start process on $(station)")

    # open file io at each process
    fi = jldopen(InputDict["RawData_path"], "r")

    bt_kurtosis = 0.0
    bt_stalta = 0.0

    S = SeisData()

    for fikey in keys(fi["Waveforms/$(station)"])

        #read SeisChannel data
        S1_raw = fi["Waveforms/$(station)/$(fikey)"]
        S1 = deepcopy(S1_raw)
        # skip if it does not exist
        if isempty(S1.x) || all(isnan.(S1.x))
            continue
        end

        #NOTE:Transient error on combination between download margin and resampling#
        # S.t may contain index 0 in t[:, 1]. In that case, this causes the error on later part, so discard that data.
        if any(S1.t[:, 1] .== 0)
            @warn("zero index found. discard it.")
            continue
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

        # Append raw trace, kurtosis and stalta time traces to S1
        if InputDict["Append_alltraces"]
            S1.misc["raw_trace"] = S1_raw.x
            S1.misc["kurtsis_trace"] = S1.misc["kurtosis"]
            S1.misc["stalta_trace"] = S1.misc["stalta"]
            delete!(S1.misc, "kurtosis")
            delete!(S1.misc, "stalta")
        else
            # don't keep those traces to reduce data size
            delete!(S1.misc, "kurtosis")
            delete!(S1.misc, "stalta")
            delete!(S1.misc, "noisesignal")
        end
        push!(S, S1)
    end


    temppath = joinpath(InputDict["tmpdir"], station * ".jld2")

    jldopen(temppath, "w") do fo
        fo["S"] = S
    end

    JLD2.close(fi)

    return (bt_kurtosis, bt_stalta)
end
