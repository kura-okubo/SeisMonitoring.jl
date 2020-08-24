"""
    scan_stations(finame::String))

scan stations saved in SeisMointoring.jl jld2 file.
"""
function scan_stations(finame::String)

    # fi = jldopen(finame, "r")
    # StationDict = OrderedDict{String, Array{String, 1}}()
    # for key in keys(fi["Waveforms"])
	#
    #     if !haskey(StationDict, key)
    #         StationDict[key] = Array{String, 1}(undef, 0)
    #     end
	#
    #     filenames = keys(fi[joinpath("Waveforms", key)])
	#
    #     for fname in filenames
    #         netstachan = split(fname, "__")[1]
    #         if  netstachan ∉ StationDict[key]
    #             push!(StationDict[key], netstachan)
    #         end
    #     end
    # end
	paths = readdir(finame)
	StationDict = OrderedDict{String, Array{String, 1}}()
	all_stations = String[]

	for path in paths
		netstachan = split(path, "__")[1]
		key = join(split(netstachan, ".")[1:2], ".")
	    if !haskey(StationDict, key)
	        StationDict[key] = Array{String, 1}(undef, 0)
	    end
		if  netstachan ∉ StationDict[key]
			push!(StationDict[key], netstachan)
			push!(all_stations, netstachan)
		end
	end

    return StationDict, all_stations
end

# finame = "RawData.jld2"
# StationDict = scan_stations(finame)
#

#
# """
#     get_stationpairs(StationDict::Dict, cc_normalization::String="cross-correlation", pairs_option::Array{String, 1}=["all"])
#
# get stations pairs to compute cross-correlations.
#
# # Argument
# - `StationDict::`: Station dictionary made by `scan_stations()`
# - `cc_normalization::String`: cross-correlation normalization method ("none", "coherence", "deconvolution")
# - `pairs_option`: combination of components (e.g. ["XX", "YY", "ZZ"] or ["all"])
#
# # Return
# - `StationPairDict::Dict` Station pair dictionary containing all station and channel pairs to be computed.
# """
# function get_stationpairs(StationDict::OrderedDict, cc_normalization::String="cross-correlation",
# 	pairs_option::Array{SubString{String},1}=["all"], chanpair_type::Array{SubString{String},1}=["all"])
#
#     stations = collect(keys(StationDict))
#     Nstation = length(stations)
#
# 	# StationPairDict = OrderedDict{String, Array{String, 1}}()
# 	StationPairs = String[]
#
# 	netstachan1_filtered = String[]
# 	netstachan2_filtered = String[]
#
#     for i = 1:Nstation
#
# 		if lowercase(cc_normalization) == "deconvolution"
# 			# compute all station pairs as deconvolution process is asymmetric.
# 			jinit = 1
#
#         elseif lowercase(cc_normalization) == "none" || lowercase(cc_normalization) == "coherence"
#             # pairwise only upper triangle of pairs
#             jinit = i
#         else
#             @error("unknown cc_normalization method: $(cc_normalization)")
#         end
#
#         for j = jinit:Nstation
#
#             sta1 = stations[i]
#             sta2 = stations[j]
#
#             # pairkey = join([sta1, sta2], "-")
#
#             # if !haskey(StationPairDict, pairkey)
# 			#
#             #     StationPairDict[pairkey] = Array{String, 1}(undef, 0)
#             # end
#
#             # evaluate components
#
#             netstachan1_all = StationDict[sta1]
#             netstachan2_all = StationDict[sta2]
#
#             for netstachan1 in netstachan1_all
#                 for netstachan2 in netstachan2_all
#
#                     #netstachan1 and netstachan2 are e.g. "BP.CCRB..BP1" and "BP.EADB..BP2"
#                     #parse components and push it into pairdict
#                     paircomp= netstachan1[end]*netstachan2[end]
# 					ct = get_chanpairtype(string.([netstachan1, netstachan2]))
#                     if (paircomp ∈ pairs_option || "all" ∈ pairs_option) && (ct ∈ chanpair_type || "all" ∈ chanpair_type)
#                         # pairs_option can be either "XX, YY, ZZ, XY..." or "all"
#                         # this pair is added to StationPairDict
# 						channelkey = join([netstachan1, netstachan2], "-")
# 						# push!(StationPairDict[pairkey], join([netstachan1, netstachan2], "-"))
# 						push!(StationPairs, channelkey)
# 						push!(netstachan1_filtered, netstachan1)
# 						push!(netstachan2_filtered, netstachan2)
#                     end
#                 end
#             end
#         end
#     end
# 	# return StationPairs
# 	return netstachan1_filtered, netstachan2_filtered, StationPairs
# end
#
# cc_normalization = "deconvolution" #"cross-correlation"
# pairs_option = String["all"]
# get_stationpairs(StationDict, cc_normalization, pairs_option)
#


"""
    get_stationpairs_chunk(all_stations_chunk::Array{String,1}, cc_normalization::String="cross-correlation", pairs_option::Array{String, 1}=["all"])

get stations pairs to compute cross-correlations.

# Argument
- `StationDict::`: Station dictionary made by `scan_stations()`
- `cc_normalization::String`: cross-correlation normalization method ("none", "coherence", "deconvolution")
- `pairs_option`: combination of components (e.g. ["XX", "YY", "ZZ"] or ["all"])

# Return
- `StationPairDict::Dict` Station pair dictionary containing all station and channel pairs to be computed.
"""
function get_stationpairs_chunk(all_stations_chunk::Array{String,1}, cc_normalization::String="cross-correlation",
	pairs_option::Array{SubString{String},1}=["all"], chanpair_type::Array{SubString{String},1}=["all"])

	#1. make StationDict
	StationDict = OrderedDict{String, Array{String, 1}}()
	all_stations = String[]

	for netstachan in all_stations_chunk
		key = join(split(netstachan, ".")[1:2], ".")
		if !haskey(StationDict, key)
			StationDict[key] = Array{String, 1}(undef, 0)
		end
		if  netstachan ∉ StationDict[key]
			push!(StationDict[key], netstachan)
		end
	end

    stations = collect(keys(StationDict))
    Nstation = length(stations)

	# StationPairDict = OrderedDict{String, Array{String, 1}}()
	StationPairs = String[]

	netstachan1_filtered = String[]
	netstachan2_filtered = String[]

    for i = 1:Nstation

		if lowercase(cc_normalization) == "deconvolution"
			# compute all station pairs as deconvolution process is asymmetric.
			jinit = 1

        elseif lowercase(cc_normalization) == "none" || lowercase(cc_normalization) == "coherence"
            # pairwise only upper triangle of pairs
            jinit = i
        else
            @error("unknown cc_normalization method: $(cc_normalization)")
        end

        for j = jinit:Nstation

            sta1 = stations[i]
            sta2 = stations[j]

            # evaluate components

            netstachan1_all = StationDict[sta1]
            netstachan2_all = StationDict[sta2]

            for netstachan1 in netstachan1_all
                for netstachan2 in netstachan2_all

                    #netstachan1 and netstachan2 are e.g. "BP.CCRB..BP1" and "BP.EADB..BP2"
                    #parse components and push it into pairdict
                    paircomp= netstachan1[end]*netstachan2[end]
					ct = get_chanpairtype(string.([netstachan1, netstachan2]))
                    if (paircomp ∈ pairs_option || "all" ∈ pairs_option) && (ct ∈ chanpair_type || "all" ∈ chanpair_type)
                        # pairs_option can be either "XX, YY, ZZ, XY..." or "all"
                        # this pair is added to StationPairDict
						channelkey = join([netstachan1, netstachan2], "-")
						push!(StationPairs, channelkey)
						push!(netstachan1_filtered, netstachan1)
						push!(netstachan2_filtered, netstachan2)
                    end
                end
            end
        end
    end
	# return StationPairs
	return netstachan1_filtered, netstachan2_filtered, StationPairs
end


function get_cc_time_windows(cc_time_unit::Int64, fs::Float64, starttime::DateTime, endtime::DateTime)

    su = d2u(starttime)
    eu = d2u(endtime)
    starts = Array(range(su,stop=eu-cc_time_unit,step=cc_time_unit))
	# ends = starts .+ cc_time_unit .- 1. / fs
	# DEBUG: overlap end[i] and start[i+1]
	ends = starts .+ cc_time_unit

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

    all_stationchan = String[]
    for stationpair in StationPairList
        sta1, sta2 = split(stationpair, "-")
        sta1 ∉ all_stationchan && push!(all_stationchan, sta1)
        sta2 ∉ all_stationchan && push!(all_stationchan, sta2)
    end
    return all_stationchan
end

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
"""
    cc_medianmute(C::CorrData, cc_medianmute_α::Float64 = 10.0)

Mute cross-correlation function whose maximum amplitude is more than
`cc_medianmute_α*median(cross-correlation functions)`
"""
function cc_medianmute!(C::CorrData, cc_medianmute_α::Float64 = 10.0)
	C.corr, inds = cc_medianmute(C.corr, cc_medianmute_α)
	C.t = remove_medianmute(C, inds)
	return nothing
end


function cc_medianmute(A::AbstractArray, cc_medianmute_α::Float64 = 10.0)

    #1. compute median of maximum amplitude of all corrs
    T, N = size(A)

    cc_maxamp = vec(maximum(abs.(A), dims=1))
    cc_medianmax = median(cc_maxamp)
    inds = findall(x-> x <= cc_medianmute_α*cc_medianmax,cc_maxamp)

    #NOTE: you cannot unbind entire array, so remove_nanandzerocol! is not used here.
    return A[:, inds], inds
end

remove_medianmute(C::CorrData, inds) = (return C.t[inds])


"""
    get_chunk_fi_stationdict(rawdata_path_all::Array{String, 1}, st::DateTime, et::DateTime)
return stationdict to copy data from shared storage to local tmp directory.
"""
function get_chunk_fi_stationdict(rawdata_path_all::Array{String, 1}, st::DateTime, et::DateTime)

	chunk_fi_stationdict = Dict{String,Array{String,1}}()
	for rawdata_path in rawdata_path_all
		finame = splitdir(rawdata_path)[2]
		netstachan, fist_str, fiet_str, _ = split(finame, "__")
		fist, fiet = DateTime.([fist_str, fiet_str])
		if (st<=fist&&fiet<=et) || (fist<st&&fiet>st) || (fist<et&&fiet>et)
			#add this file into list
			!haskey(chunk_fi_stationdict, netstachan) && (chunk_fi_stationdict[netstachan] = String[])
			push!(chunk_fi_stationdict[netstachan], finame)
		end
	end
	return chunk_fi_stationdict
end

# """
#     precompile_map_compute_fft(all_stations::Array{String,1}, InputDict::Dict)
# Call `map_compute_fft()` in all workers to precompile.
#
# # NOTE
#
# The number of required workers for `compute_fft` is generally smaller than that of `compute_cc`.
# For example, `N_station` requires `N_station` workers for fft and `N_station*(N_station-1)` for cross-correlation pair.
# When using large number of cores as workers to parallelize `compute_cc` processes,
# the number of workers is too much for `compute_fft` process.
#
# This causes an issue such that the asyncronous call of `compute_fft` randomly increases precompile time.
#
# e.g. Given 100 cores as workers, and parallelize 10 stations for fft:
#
# - first chunk: worker 1,2 ... 10 precompiled
# - second chunk: worker 1,2, ..., 9, 11, taking precompile time for worker 11
# - third chunk: worker 1,2, ..., 9, 12, taking precompile time for worker 12
# - fourth chunk: worker 1,2, ..., 9, 13, taking precompile time for worker 13
#
# This function calls `compute_fft` in the beginning to avoid the redundant precompiling.
#
# """
# function precompile_map_compute_fft(rawdata_path_all::Array{String,1}, all_stations::Array{String,1}, InputDict::OrderedDict)
#
# 	println("---start precompile map_compute_fft---")
#
# 	# # processing short time chunk
# 	# InputDict["starts_chunk"] = InputDict["starts"][1:3]
# 	# InputDict["ends_chunk"] = InputDict["ends"][1:3]
# 	#
# 	# # make pseudo station name list to allocate compute_fft to all workers
# 	# k = ceil(Int, nworkers()/length(all_stations))
# 	# precompile_stationlist = repeat(all_stations, outer=k)
# 	# println("debug: length precompile_stationlist = $(length(precompile_stationlist[1:nworkers()])).")
# 	#
# 	# if InputDict["use_local_tmpdir"]
# 	# 	# make chunk_fi_stationdict to copy data to local tmp directory
# 	# 	InputDict["chunk_fi_stationdict"] = get_chunk_fi_stationdict(rawdata_path_all,
# 	# 											u2d(InputDict["starts_chunk"][1]), u2d(InputDict["ends_chunk"][end]))
# 	# end
# 	#
# 	# t_precompilefft = @elapsed pmap(x -> map_compute_fft(x, InputDict, precompile_mode=true), precompile_stationlist[1:nworkers()])
#
# 	println("---precompile map_compute_fft done with $(t_precompilefft)[s]---")
#
# end
# key = "BP.CCRB-BP.EADB"
# all_stationchannels = get_stationchanname(StationPairDict[key])

#
# """
# 	remove_nancol(A::AbstractArray)
#
# Remove column (i.e. trace) which has NaN.
# """
# function remove_nanandzerocol(A::AbstractArray)
#
# 	N = size(A, 2)
# 	nancol = ones(Int64, N)
# 	for i = 1:N
# 		if any(isnan.(A[:, i])) || all(iszero, A[:,i])
# 			# this has NaN in its column
# 			nancol[i] = 0
# 		end
# 	end
#
# 	nancol=convert.(Bool, nancol)
#
# 	#NOTE: you cannot unbind entire array, so remove_nanandzerocol! is not used here.
# 	return A[:, nancol], nancol
#
# end
#
# remove_nanandzerocol_t(C, nancol) = (return C.t[nancol])
#
# """
# 	remove_nanandzerocol!(C::CorrData)
#
# Remove any nan and all zero column from CorrData.
# Modify C.t as well.
# """
#
# function remove_nanandzerocol!(C::CorrData)
# 	C.corr, nancol = remove_nanandzerocol(C.corr)
# 	C.t = remove_nanandzerocol_t(C, nancol)
# end
