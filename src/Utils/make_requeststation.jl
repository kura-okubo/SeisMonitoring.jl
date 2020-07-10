using DataFrames, Dates, JLD2, CSV

"""
    make_requeststation_fromIRISgmap(fopath::String)

make request stations from IRIS Gmap text output with Julia DataFrames format.

Request station dataframe format:
DataFrame(network=String[],
                           station=String[],
                           location=String[],
                           channel=String[],
                           latitude=String[],
                           longitude=String[],
                           elevation=String[])
for each earthquake data center.
We call it as 'request_src_chunk'.

Then, fo.jld2 contains each src chunk such as:
./default_requeststations.jld2
 └─ NCEDC
 └─ IRIS

See available src chunk name in SeisIO.seis_www
"""
function make_requeststation_fromIRISgmap(fipath::String; locchan::Dict=Dict(), fodir::String="", foname::String="request_stations")

    # default fo path is same with fipath directory
    isempty(fodir) && (fodir = splitdir(fipath)[1])
    fopath = joinpath(fodir, foname)
    foname[end-4:end] != ".jld2" && (fopath *= ".jld2")

    # read input file
    lines = readlines(fipath)

    # parse algorithm:
    # 1. search until encounting #DATACENTER
    # 2. append stations to the chunc

    RequestStations = Dict()
    i = 1
    while i <= length(lines)
        if occursin("#DATACENTER", lines[i])
            #start appending this srcchunk
            src = split(split(lines[i], "=")[2], ",")[1]
            # #manipulate IRIS variation
            # occursin("IRIS", src) && (src = "IRIS")
            srcchunk = DataFrame(network=String[],
                                   station=String[],
                                   location=String[],
                                   channel=String[],
                                   latitude=String[],
                                   longitude=String[],
                                   elevation=String[],
                                   sitename=String[],
                                   starttime=DateTime[],
                                   endtime=DateTime[])

            i += 1
            while true
                # append request stations
                lc = split(lines[i],"|")
                net=lc[1]
                sta=lc[2]
                lat=lc[3]
                lon=lc[4]
                ele=lc[5]
                sitename=lc[6]
                starttime=DateTime(lc[7])
                endtime=DateTime(lc[8])

                # choose request channel if locchan dict contains the network; otherwise request all (e.g. *, *)
                if net ∈ keys(locchan)
                    reqs = locchan[net]
                    for req in reqs
                        loc, chan = string.(req)

                        dftemp = DataFrame(network=net, station=sta,location=loc,channel=chan,
                        latitude=lat, longitude=lon,elevation=ele, sitename=sitename, starttime=starttime,endtime=endtime)
                        append!(srcchunk,dftemp)
                    end

                else
                    dftemp = DataFrame(network=net, station=sta,location="*",channel="*",
                    latitude=lat, longitude=lon,elevation=ele, sitename=sitename, starttime=starttime,endtime=endtime)
                    append!(srcchunk,dftemp)
                end

                i += 1
                if (i > length(lines) || isempty(lines[i]))
                    # push srcchunk to master Dictionary
                    RequestStations[src] = srcchunk

                    # output csv
                    CSV.write(fopath[1:end-5]*"_$(src).csv", srcchunk)
                    # next src chunk or end of file
                    break
                end
            end

        else
            i += 1
            continue
        end

    end

    # output jld2 file for SeisMonitoring
    jldopen(fopath, "w") do f
        for key in keys(RequestStations)
            f[key] = RequestStations[key]
        end
    end

end


# fipath = "/Users/kurama/Dropbox/BP_monitoring/dev_make_requeststation/localtestBP2004_tmp_INPUT/BP_gmap-stations.txt"
# # fipath = "/Users/kurama/Dropbox/BP_monitoring/dev_make_requeststation/testgmap-stations-3.txt"
#
# locchan = Dict(
#             "BP" => [(40, "L*I")]
#             )
#
# make_requeststation_fromIRISgmap(fipath, locchan=locchan, foname="BP_temperature.jld2")
