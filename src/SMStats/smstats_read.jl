using SeisIO, SeisNoise, Dates, JLD2, DataFrames, CSV, Distributed
"""
    smstats_read()

Get the statistics of change of velocity and attenuation history, and save them to CSV file.


"""
function smstats_read(shorttimestackdir::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    foname::String = "monitoring_stats.csv")

    !ispath(fodir) && mkpath(fodir)

    println("-------Read dvv and dQc--------")

    paths = SeisIO.ls(shorttimestackdir)

    df_mapped = pmap(x -> map_get_monitoringdf(x, starttime, endtime), paths)
    #df_mapped = pmap(x -> map_get_monitoringdf(x, starttime, endtime), paths[1:10])

    # append the data frame

    df_all = DataFrame(date=DateTime[], stationpair=String[], networks=String[], components=String[], freqband=String[],
                        cc_dvv=Float64[], dvv=Float64[], dist_dualstretch=Float64[], dQcinv = Float64[], dAA=Float64[])

    for df = df_mapped
        !isnothing(df) && append!(df_all, df)
    end

    csvname = joinpath(fodir, foname)
    CSV.write(csvname, df_all)

    @info("$(splitdir(csvname)[2]) is successfully saved. For plotting, please use smplot_dvv(\"$(csvname)\", kwargs...) or smplot_dQc($(csvname), kwargs...)")

    return nothing
end

function map_get_monitoringdf(path::String, starttime::DateTime, endtime::DateTime)

    fi = try
        jldopen(path, "r")
    catch
        # file cannot be opened
        return nothing
    end

    # loop keys


    df_path = DataFrame(date=DateTime[], stationpair=String[], networks=String[], components=String[], freqband=String[],
                        cc_dvv=Float64[], dvv=Float64[], dist_dualstretch=Float64[], dQcinv = Float64[], dAA=Float64[])

    for stachankey in keys(fi)
        println("start reading $(stachankey)")
        for timekey in keys(fi[stachankey])
            for freqkey in keys(fi[joinpath(stachankey, timekey)])
                #1. parse metadata
                station1, station2 = split(stachankey, "-")
                networks = join([split(station1, ".")[1], split(station2, ".")[1]], "-")
                components = station1[end]*station2[end]

                st, et = DateTime.(split(timekey, "--"))
                if !(starttime >= et || endtime <= st)
                # this file has overlap with the target timewindow [starttime, endtime]
                    corrkey = joinpath(stachankey, timekey, freqkey)
                    C = fi[corrkey]
                    if !haskey(C.misc, "dvv") && !haskey(C.misc, "dQcinv")
                        @warn("$(corrkey) does not have dvv nor dQcinv. Please check measurement_method in the input file.")
                        continue;
                    end

                    # append dataframe
                    df = DataFrame(date = C.misc["stack_centraltime"], stationpair=stachankey,
                                    networks=networks, components=components, freqband=freqkey)
                    # append dvv
                    if haskey(C.misc, "dvv")
                        df_dvv = DataFrame(date = C.misc["stack_centraltime"], cc_dvv = float(C.misc["cc_dvv"]), dvv = float(C.misc["dvv"]))
                        df = leftjoin(df, df_dvv, on=:date)
                    end

                    # append dQc and dAA
                    if haskey(C.misc, "dQcinv")
                        df_dQc = DataFrame(date = C.misc["stack_centraltime"], dist_dualstretch = float(C.misc["dist_dualstretch"]),
                         dQcinv = float(C.misc["dQcinv"]), dAA = float(C.misc["dAA"]))
                        df = leftjoin(df, df_dQc, on=:date)
                    end

                    append!(df_path, df)
                end
            end
        end
    end
    close(fi)

    return df_path
end
