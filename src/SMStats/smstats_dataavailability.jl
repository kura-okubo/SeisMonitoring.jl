using SeisIO, CSV
"""
    smstats_dataavailability(fidir::String,foname::String,
                            starttime::DateTime,
                            endtime::DateTime;
                            network::Union{String, AbstractArray}=["all"]
                            )
Compute noise data fraction as data abvailability.
This function is extended from smplot_noiseavailability.jl

# Argument
- `fidir::String`: absolute/relative path to seismicdata with seisio format.
- `foname::String`: figure and csv output directory.
- `starttime::DateTime`: starttime to be plotted.
- `endtime::DateTime`: endtime to be plotted.
- `network::Union{String, AbstractArray}=["all"]`: Array of network to be plotted. (e.g. ["all"], or ["BP, "NC"])
"""
function smstats_dataavailability(fidir::String,foname::String,
    starttime::DateTime,
    endtime::DateTime;
    network::Union{String, AbstractArray}=["all"])

    typeof(network) == String && (network = [network])

    PlotDict=Dict(
    "fidir" => fidir,
    "foname" => foname,
    "starttime" => starttime,
    "endtime" => endtime,
    "network" => network
    )

    println("***************************************")
    println("Plot Noise Availability")
    println("starttime    = $(starttime)")
    println("endtime      = $(endtime)")
    println("network      = $(network)")
    println("***************************************\n")

	stations = String[]
	for (root, dirs, files) in walkdir(PlotDict["fidir"])
       for file in files
		   fi = joinpath(root, file)
		   (split(fi, ".")[end] == "seisio") && push!(stations, fi)# filter if it is .seisio
       end
    end

    # filter the network
    "all" ∉ network && filter!(x -> split(splitdir(x)[2], ".")[1] ∈ network, stations)

    # parallelize with keys in Rawdata.jld2 i.e. stations
    println("-------START Getting Noise Data Fraction--------")

    df_mapped = pmap(x -> map_getnoisedatafraction_dataavail(x, PlotDict), stations)

    df_all = DataFrame(station=String[], date = String[], data_fraction=Float64[], removal_fraction=Float64[])

    for df in df_mapped
        !isnothing(df) && append!(df_all, df);
    end

    # output csvfile

    CSV.write(foname*".csv", df_all)

    println("smstats_dataavailability is successfully done.")

    return nothing
end

function map_getnoisedatafraction_dataavail(stationpath::String, PlotDict::Dict)

    fidir, fikey = splitdir(stationpath)

	# process only .seisio file
	if split(fikey, ".")[end] != "seisio"
		return nothing
	end

    # println("start process on $(fikey)")

    df_station = DataFrame(station=String[], date = String[], data_fraction=Float64[], removal_fraction=Float64[])

    # fi = jldopen(PlotDict["filename"], "r")

    # for fikey in keys(fi["Waveforms/$(station)"])
    file_st, file_et = DateTime.(split(fikey, "__")[2:3])

    if !(PlotDict["starttime"] >= file_et || PlotDict["endtime"] <= file_st)
        # this file has overlap with the target timewindow [starttime, endtime]
		# S1 = fi["Waveforms/$(station)/$(fikey)"]
		S1 = rseis(stationpath)[1]
        centraltime = string(u2d((d2u(file_st) + d2u(file_et))/2))
        data_fraction = S1.misc["data_fraction"]

		if haskey(S1.misc, "removal_fraction")
			# this data is for remove seisEQ
			removal_fraction = S1.misc["removal_fraction"]
		else
			removal_fraction = 0.0 # raw data (no removal)
		end

        push!(df_station, (S1.id, centraltime, data_fraction, removal_fraction))
    end

    return df_station

end
