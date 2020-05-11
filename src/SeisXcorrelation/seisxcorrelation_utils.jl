"""
    scan_stations(finame::String))

scan stations saved in SeisMointoring.jl jld2 file.
"""
function scan_stations(finame::String)

    fi = jldopen(finame, "r")
    StationDict = OrderedDict{String, Array{String, 1}}()
    for key in keys(fi["Waveforms"])

        if !haskey(StationDict, key)
            StationDict[key] = Array{String, 1}(undef, 0)
        end

        filenames = keys(fi[joinpath("Waveforms", key)])

        for fname in filenames
            netstachan = split(fname, "__")[1]
            if  netstachan ∉ StationDict[key]
                push!(StationDict[key], netstachan)
            end
        end
    end
    return StationDict
end

# finame = "RawData.jld2"
# StationDict = scan_stations(finame)
#


"""
    get_stationpairs(StationDict::Dict, cc_method::String="cross-correlation", pairs_option::Array{String, 1}=["all"])

get stations pairs to compute cross-correlations.

# Argument
- `StationDict::`: Station dictionary made by `scan_stations()`
- `cc_method::String`: cross-correlation method ("cross-correlation", "coherence", "deconvolution")
- `pairs_option`: combination of components (e.g. ["XX", "YY", "ZZ"] or ["all"])

# Return
- `StationPairDict::Dict` Station pair dictionary containing all station and channel pairs to be computed.
"""
function get_stationpairs(StationDict::OrderedDict, cc_method::String="cross-correlation", pairs_option::Array{SubString{String},1}=["all"])

    stations = collect(keys(StationDict))
    Nstation = length(stations)

    StationPairDict = OrderedDict{String, Array{String, 1}}()

    for i = 1:Nstation

        if lowercase(cc_method) == "cross-correlation" || lowercase(cc_method) == "coherence"
            # pairwise only upper triangle of pairs
            jinit = i
        elseif lowercase(cc_method) == "deconvolution"
            # compute all station pairs as deconvolution process is asymmetric.
            jinit = 1
        else
            @error("unknown cc-method: $(cc_method)")
        end

        for j = jinit:Nstation

            sta1 = stations[i]
            sta2 = stations[j]

            pairkey = join([sta1, sta2], "-")

            if !haskey(StationPairDict, pairkey)

                StationPairDict[pairkey] = Array{String, 1}(undef, 0)
            end

            # evaluate components

            netstachan1_all = StationDict[sta1]
            netstachan2_all = StationDict[sta2]

            for netstachan1 in netstachan1_all
                for netstachan2 in netstachan2_all

                    #netstachan1 and netstachan2 are e.g. "BP.CCRB..BP1" and "BP.EADB..BP2"
                    #parse components and push it into pairdict
                    paircomp= netstachan1[end]*netstachan2[end]

                    if paircomp ∈ pairs_option || "all" ∈ pairs_option
                        # pairs_option can be either "XX, YY, ZZ, XY..." or "all"
                        # this pair is added to StationPairDict
                        push!(StationPairDict[pairkey], join([netstachan1, netstachan2], "-"))
                    end
                end
            end
        end
    end
    return StationPairDict
end
#
# cc_method = "deconvolution" #"cross-correlation"
# pairs_option = String["all"]
# get_stationpairs(StationDict, cc_method, pairs_option)
#

function get_cc_time_windows(cc_time_unit::Int64, fs::Float64, starttime::DateTime, endtime::DateTime)

    su = d2u(starttime)
    eu = d2u(endtime)
    starts = Array(range(su,stop=eu,step=cc_time_unit))
    ends = starts .+ cc_time_unit .- 1. / fs

    return starts, ends
end

#
# cc_time_unit=86400 #[sec]
# fs=20.0 #[sec]
# starttime = DateTime("2004-04-01T00:00:00")
# endtime = DateTime("2004-04-03T00:00:00")
# starts, ends = get_cc_time_windows(cc_time_unit, fs, starttime, endtime)


"""
    get_stationchanname(StationPairList::)
get all station channel name used in the station pair list.

# Argument
- `StationPairList::Array{String, 1}` List of stastion pairs from get_stastionpairs.jl

e.g.    get_stationchanname(StationPairDict[key])

# Return

- `StationChan::Array{String,1}` List of stations used in StationPairList.
"""
function get_stationchanname(StationPairList::Array{String,1})

    all_stationchan = []
    for stationpair in StationPairList
        sta1, sta2 = split(stationpair, "-")
        sta1 ∉ all_stationchan && push!(all_stationchan, sta1)
        sta2 ∉ all_stationchan && push!(all_stationchan, sta2)
    end
    return all_stationchan
end


# key = "BP.CCRB-BP.EADB"
# all_stationchannels = get_stationchanname(StationPairDict[key])


"""
	remove_nancol(A::AbstractArray)

Remove column (i.e. trace) which has NaN.
"""
function remove_nanandzerocol(A::AbstractArray)

	N = size(A, 2)
	nancol = ones(Int64, N)
	for i = 1:N
		if any(isnan.(A[:, i])) || all(iszero, A[:,i])
			# this has NaN in its column
			nancol[i] = 0
		end
	end

	nancol=convert.(Bool, nancol)

	#NOTE: you cannot unbind entire array, so remove_nanandzerocol! is not used here.
	return A[:, nancol], nancol

end

remove_nanandzerocol_t(C, nancol) = (return C.t[nancol])

"""
	remove_nanandzerocol!(C::CorrData)

Remove any nan and all zero column from CorrData.
Modify C.t as well.
"""

function remove_nanandzerocol!(C::CorrData)
	C.corr, nancol = remove_nanandzerocol(C.corr)
	C.t = remove_nanandzerocol_t(C, nancol)
end
