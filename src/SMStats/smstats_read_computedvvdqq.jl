using SeisIO, SeisNoise, Dates, JLD2, DataFrames, CSV, Distributed
"""
    smstats_read_computedvvdqq()

Get the statistics of change of velocity and attenuation history, and save them to CSV file.


"""
function smstats_read_computedvvdqq(shorttimestackdir::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    foname::String = "monitoring_stats.csv")

    !ispath(fodir) && mkpath(fodir)

    println("-------Read dvv and dQc--------")

    paths = SeisIO.ls(shorttimestackdir)

    df_mapped = pmap(x -> map_get_monitoringdf_comnputedvvdqq(x, starttime, endtime), paths)
    #df_mapped = pmap(x -> map_get_monitoringdf(x, starttime, endtime), paths[1:10])

    # append the data frame

    df_all = DataFrame(date=DateTime[], stationpair=String[], networks=String[], components=String[], freqband=String[],
                        cc_dvv=Float64[], dvv=Float64[],
                        dqq_pos=Float64[], dqq_neg=Float64[], dqq_avg = Float64[],
                        dss_pos=Float64[], dss_neg=Float64[], dss_avg = Float64[],
                        Qcinv_pos_ref=Float64[],Qcinv_neg_ref=Float64[],Qcinv_pos_cur=Float64[],Qcinv_neg_cur=Float64[])

    for df = df_mapped
        !isnothing(df) && append!(df_all, df)
    end

    csvname = joinpath(fodir, foname)
    CSV.write(csvname, df_all)

    @info("$(splitdir(csvname)[2]) is successfully saved. For plotting, please use smplot_dvv(\"$(csvname)\", kwargs...) or smplot_dQc($(csvname), kwargs...)")

    return nothing
end

function map_get_monitoringdf_comnputedvvdqq(path::String, starttime::DateTime, endtime::DateTime)

    fi = try
        jldopen(path, "r")
    catch
        # file cannot be opened
        return nothing
    end

    # loop keys
    df_path = DataFrame(date=DateTime[], stationpair=String[], networks=String[], components=String[], freqband=String[],
                        cc_dvv=Float64[], dvv=Float64[],
                        dqq_pos=Float64[], dqq_neg=Float64[], dqq_avg = Float64[],
                        dss_pos=Float64[], dss_neg=Float64[], dss_avg = Float64[],
                        Qcinv_pos_ref=Float64[],Qcinv_neg_ref=Float64[],Qcinv_pos_cur=Float64[],Qcinv_neg_cur=Float64[])

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
                    if !haskey(C.misc, "dvv") && !haskey(C.misc, "dqq_avg")
                        @warn("$(corrkey) does not have dvv nor dqq_avg. Please check measurement_method in the input file.")
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

                    # append dqq and dss
                    if haskey(C.misc, "dqq_avg")
                        df_dqq = DataFrame(date = C.misc["stack_centraltime"],
                        dqq_pos= float(C.misc["dqq_pos"]), dqq_neg= float(C.misc["dqq_neg"]), dqq_avg = float(C.misc["dqq_avg"]))
                        df = leftjoin(df, df_dqq, on=:date)

                        df_dss = DataFrame(date = C.misc["stack_centraltime"],
                        dss_pos= float(C.misc["dss_pos"]), dss_neg= float(C.misc["dss_neg"]), dss_avg = float(C.misc["dss_avg"]))
                        df = leftjoin(df, df_dss, on=:date)
                    end

                    # append Qc
                    if haskey(C.misc, "Qcinv_pos_cur")
                        df_Qcinv = DataFrame(date = C.misc["stack_centraltime"],
                            Qcinv_pos_ref= float(C.misc["Qcinv_pos_ref"]), Qcinv_neg_ref= float(C.misc["Qcinv_neg_ref"]),
                            Qcinv_pos_cur= float(C.misc["Qcinv_pos_cur"]), Qcinv_neg_cur= float(C.misc["Qcinv_neg_cur"]),)
                        df = leftjoin(df, df_Qcinv, on=:date)

                    end


                    append!(df_path, df)
                end
            end
        end
    end
    close(fi)

    return df_path
end
