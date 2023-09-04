using SeisIO, SeisNoise, Dates, JLD2, DataFrames, CSV, Distributed
"""
    smstats_read_stretching()

Get the statistics of change of velocity history, and save them to CSV file.


"""
function smstats_read_stretching(shorttimestackdir::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    foname::String = "monitoring_stats.csv")

    !ispath(fodir) && mkpath(fodir)

    println("-------Read dvv--------")

    paths = SeisIO.ls(shorttimestackdir)

    df_mapped = pmap(x -> map_get_monitoringdf_stretching(x, starttime, endtime), paths)
    #df_mapped = pmap(x -> map_get_monitoringdf(x, starttime, endtime), paths[1:10])

    # append the data frame
    df_all = DataFrame(date=DateTime[], stationpair=String[], networks=String[], components=String[], freqband=String[],
                        dvv_ts=Float64[], cc_ts=Float64[], err_ts=Float64[])

    for df = df_mapped
        !isnothing(df) && append!(df_all, df)
    end

    csvname = joinpath(fodir, foname)
    CSV.write(csvname, df_all)

    @info("$(splitdir(csvname)[2]) is successfully saved.)")

    return nothing
end

function map_get_monitoringdf_stretching(path::String, starttime::DateTime, endtime::DateTime)

    fi = try
        jldopen(path, "r")
    catch
        # file cannot be opened
        return nothing
    end

    # loop keys
    df_path = DataFrame(date=DateTime[], stationpair=String[], networks=String[], components=String[], freqband=String[],
                        dvv_ts=Float64[], cc_ts=Float64[], err_ts=Float64[])

    # for stachankey in keys(fi)
    stachankey = split(splitdir(path)[2][1:end-5], "_")[2]
    println("start reading $(stachankey)")
    for timekey in keys(fi)
        for freqkey in keys(fi[timekey])
            #1. parse metadata
            # station1, station2 = split(stachankey, "-")
            # networks = join([split(station1, ".")[1], split(station2, ".")[1]], "-")
            # components = station1[end]*station2[end]
            station1, station2, components = split(stachankey, "-")
            networks = join([split(station1, ".")[1], split(station2, ".")[1]], "-")

            st, et = DateTime.(split(timekey, "--"))

            if !(starttime >= et || endtime <= st)
            # this file has overlap with the target timewindow [starttime, endtime]
                corrkey = joinpath(timekey, freqkey)
                C = fi[corrkey]

                if !haskey(C.misc, "dvv_ts")
                    # @warn("$(corrkey) does not have dvv nor dqq_avg. Please check measurement_method in the input file.")
                    continue;
                end

                # append dataframe
                df = DataFrame(date = C.misc["stack_centraltime"], stationpair="$(station1)-$(station2)",
                                networks=networks, components=components, freqband=freqkey)
                # append dvv
                if haskey(C.misc, "dvv_ts")
                    df_dvv = DataFrame(date = C.misc["stack_centraltime"], dvv_ts = float(C.misc["dvv_ts"]), cc_ts = float(C.misc["cc_ts"]),
                                        err_ts = float(C.misc["err_ts"]))
                    df = leftjoin(df, df_dvv, on=:date)
                end

                append!(df_path, df)
            end
        end
    end
    # end
    close(fi)

    return df_path
end
