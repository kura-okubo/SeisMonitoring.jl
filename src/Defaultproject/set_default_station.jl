using DataFrames, JLD2, CSV

"""
    set_default_station(fopath::String)

make request stations in Julia DataFrames format.

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
function set_default_station(fopath::String)

    NCEDC_srcchunk = DataFrame(network=String[],
                               station=String[],
                               location=String[],
                               channel=String[],
                               latitude=String[],
                               longitude=String[],
                               elevation=String[])

    location =[""]
    channel = ["BP*"]

    module_path=splitdir(pathof(SeisMonitoring))[1]
    lines = readlines(joinpath(module_path, "Defaultproject/BP_gmap-stations.txt"))

    for line in lines
        if occursin("#", line) || isempty(line)
            #this line is commented out or empty line
            continue;
        end
        lc = split(line,"|")
        net=lc[1]
        sta=lc[2]
        lat=lc[3]
        lon=lc[4]
        ele=lc[5]

        for loc in location
            for cha in channel
                dftemp = DataFrame(network=net, station=sta,location=loc,channel=cha,
                latitude=lat, longitude=lon,elevation=ele)
                append!(NCEDC_srcchunk,dftemp)
            end
        end
    end

    RequestStations = Dict(
        "NCEDC" => NCEDC_srcchunk
    )

    CSV.write(fopath[1:end-4]*".csv", NCEDC_srcchunk)
    jldopen(fopath, "w") do f
        for key in keys(RequestStations)
            f[key] = RequestStations[key]
        end
    end
end
