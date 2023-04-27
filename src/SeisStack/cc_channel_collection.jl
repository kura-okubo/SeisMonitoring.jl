# using SeisIO, SeisNoise, JLD2, Distributed
function cc_channel_collection(ccdir::String)

    cc_outdir = joinpath(splitdir(ccdir)[1], "cc_channel_collection")
    !ispath(cc_outdir) && mkdir(cc_outdir)

    # collect station pairname
    station_pairdict = get_collect_station_pairdict(ccdir)
    # reassemble station pairs
    pmap(x -> map_collect_stationpairs(x, cc_outdir, station_pairdict), collect(keys(station_pairdict)))

    return nothing

end

function get_collect_station_pairdict(ccdir::String)

    # cc_paths = SeisIO.ls(ccdir)
    cc_paths = SeisIO.ls(ccdir*"/*") #2022.07.11 update: search jld2 files recursively
    station_pairdict = Dict{String,Array{String,1}}()
    for cc_path in cc_paths
        split(cc_path, ".")[end] != "jld2" && continue
        fname = splitdir(cc_path)[2][1:end-5]
        stachanpair, _, _ = split(fname, "__")
        stachan1, stachan2 = split(stachanpair, "-")
        net1, sta1, loc1, cha1 = split(stachan1, ".")
        net2, sta2, loc2, cha2 = split(stachan2, ".")
        nochan_stationpair = joinpath("$(net1).$(sta1)-$(net2).$(sta2)-$(cha1[end])$(cha2[end])")

        # if nochan_stationpair is not in the dict, add as new key
        if !haskey(station_pairdict, nochan_stationpair)
            station_pairdict[nochan_stationpair] = [cc_path]
        # if nochan_stationpair is in dict, but fname is not, push the key
        elseif !(fname in station_pairdict[nochan_stationpair])
            push!(station_pairdict[nochan_stationpair], cc_path)
        end
    end
    return station_pairdict
end

"""
    map_collect_stationpairs

collect corrdata into a jld2 file with respect stationpairkeys.

#Note

    - Skip the second one if two stationchannel pair have the same time stamp. (e.g. XX.XX..HHZ, XX.XX..SHZ)

"""
function map_collect_stationpairs(station_pairkey::String, cc_outdir::String, station_pairdict::Dict)

    fo = jldopen(joinpath(cc_outdir, station_pairkey*".jld2"), "w", iotype=IOStream)

    cc_paths = station_pairdict[station_pairkey] # cc jld2 files to be collected

    for cc_path in cc_paths
        fi = jldopen(cc_path, "r")
        timestamps = keys(fi)
        for timestamp in timestamps
            # check if the timestamp is in output file
            !haskey(fo, timestamp) ? JLD2.Group(fo, timestamp) : continue

            freqs = keys(fi[timestamp])
            for freq in freqs
                g1 = joinpath(timestamp, freq)
                !haskey(fo, g1) && (fo[g1] = fi[g1])
            end
        end
        close(fi)
    end
    close(fo)

    return nothing

end


# ccdir = "./cc"
# cc_channel_collection(ccdir)
