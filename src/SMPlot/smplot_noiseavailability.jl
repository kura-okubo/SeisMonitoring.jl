using VegaLite, FileIO, CSV
"""
    smplot_noiseavailability(filename::String,fodir::String,
                            starttime::DateTime,
                            endtime::DateTime;
                            network::Union{String, AbstractArray}=["all"],
                            figsize::Tuple=(1200, 800),
                            fmt::String="png"
                            )
Plot noise data fraction as data abvailability with Lasagna plot style.

# Argument
- `filename::String`: absolute/relative path to seismicdata jld2 file with SeisMonitoring.jl format.
- `fodir::String`: figure and csv output directory.
- `starttime::DateTime`: starttime to be plotted.
- `endtime::DateTime`: endtime to be plotted.
- `network::Union{String, AbstractArray}=["all"]`: Array of network to be plotted. (e.g. ["all"], or ["BP, "NC"])
- `xtimeaxisformat::String`: x label time axis format (e.g. "%Y-%b-%dT%H:%M:%S")
- `xtimeunit::String`: x label time unit to avoid too dense xlabel (e.g. "year", "month", "date")
- `xlabelrotation`: x time axis rotation
- `figsize::Tuple=(1200, 800)`: figure size to be plotted.
- `fmt::String="png"`: figure format to be plotted.
"""
function smplot_noiseavailability(filename::String,fodir::String,
    starttime::DateTime,
    endtime::DateTime;
    network::Union{String, AbstractArray}=["all"],
    xtimeaxisformat::String="%Y-%b-%dT%H:%M:%S",
    xlabelrotation::Real=0,
    xtimeunit::String = "date",
    figsize::Tuple=(1200, 800),
    fmt::String="png"
    )

    typeof(network) == String && (network = [network])

    PlotDict=Dict(
    "filename" => filename,
    "fodir" => fodir,
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

    ispath(PlotDict["filename"]) ? (fi = jldopen(PlotDict["filename"], "r")) : error("$(PlotDict["filename"]) is not found.")
    !haskey(fi, "Waveforms") && error("$(PlotDict["filename"]) does not have Waveforms group. Please check the waveform data format in JLD2.")
    stations = keys(fi["Waveforms"])
    JLD2.close(fi)

    # filter the network
    "all" ∉ network && filter!(x -> split(x, ".")[1] ∈ network, stations)

    # parallelize with keys in Rawdata.jld2 i.e. stations
    println("-------START Getting Noise Data Fraction--------")

    df_mapped = pmap(x -> map_getnoisedatafraction(x, PlotDict), stations)

    df_all = DataFrame(station=String[], date = String[], data_fraction=Float64[])

    for df in df_mapped
        append!(df_all, df);
    end

    # parse time to vegalite time format
    # NOTE: "01 Jan 2012 23:00:00" is better format for vegalite.jl
    function reformat_time(t::String)
        ymd, hms = split(t, "T")
        y, m, d = split(ymd, "-")
        return "$(d) $( Dates.monthname(parse(Int,m))) $(y) $(hms)"
    end

    df_all.date = map(x -> reformat_time(x), df_all.date)

    # println(xtimeaxisformat)

    figname = join([join(network, "_"), string(starttime), string(endtime)], "--")
    # output csvfile

    !ispath(fodir) && mkpath(fodir)

    CSV.write(joinpath(fodir, figname*".csv"), df_all)
    # plot data contents with Lasagna plot style: see VegaLite.jl document

    println("-------Plot Noise Data Availability--------")

    p = sm_vegalite_lasagnaplot(df_all, figsize=figsize,
                xtimeaxisformat=xtimeaxisformat,xlabelrotation=xlabelrotation,
                xtimeunit=xtimeunit)

    p |> FileIO.save(joinpath(fodir, figname*".$(fmt)"))

    println("smplot_noiseavailability is successfully done.")

    return nothing
end

"""
map_getnoisedatafraction(station::String, InputDict::Dict)

Return DataFrame containing stationchannel, date and noise data fraction.
"""
function map_getnoisedatafraction(station::String, PlotDict::Dict)

    println("start process on $(station)")

    df_station = DataFrame(station=String[], date = String[], data_fraction=Float64[])

    fi = jldopen(PlotDict["filename"], "r")

    for fikey in keys(fi["Waveforms/$(station)"])
        file_st, file_et = DateTime.(split(fikey, "__")[2:3])
        if !(PlotDict["starttime"] >= file_et || PlotDict["endtime"] <= file_st)
            # this file has overlap with the target timewindow [starttime, endtime]
            S1 = fi["Waveforms/$(station)/$(fikey)"]
            centraltime = string(u2d((d2u(file_st) + d2u(file_et))/2))
            data_fraction = S1.misc["data_fraction"]
            push!(df_station, (S1.id, centraltime, data_fraction))
        end
    end

    JLD2.close(fi)

    # println(df_station)

    return df_station

end


"""
sm_vegalite_lasagnaplot(df_all::DataFrame)

Plot lasagnaplot using Vegalite.jl.
Documentation can be found at:

- Vegalite.jl website (https://www.queryverse.org/VegaLite.jl/stable/)
- Offitial Vegalite webpage (https://vega.github.io/vega-lite/)
"""
function sm_vegalite_lasagnaplot(
    df_all::DataFrame;
    figsize::Tuple = (800, 600),
    xtimeaxisformat::String = "%Y-%b-%dT%H:%M:%S",
    xlabelrotation::Real = 0,
    xtimeunit::String = "date"
)

    # NOTE: So far condition of xtick span is hard to implement, so write with respect to xtimeunit.

    if lowercase(xtimeunit)=="year"
        p = (df_all |> @vlplot(
            width = figsize[1],
            height = figsize[2],
            mark = :rect,
            x = {
                field = :date,
                title = "Time",
                timeUnit = :datemonthyearhoursminutesseconds,
                axis = {
                    format = xtimeaxisformat,
                    labelAngle = xlabelrotation,
                    labelOverlap = false,
                    #===here is the custumize of x axis span==#
                    labelColor={
                        condition={
                            test={
                                field="value",
                                timeUnit=:monthdate,
                                equal={month=1,date=1}
                            },
                            value="black"
                        },
                        value=nothing
                    },
                    tickColor={
                        condition={
                            test={
                                field="value",
                                timeUnit=:monthdate,
                                equal={month=1,date=1}
                            },
                            value="black"
                        },
                        value=nothing
                    },
                    #================#
                },
            },

            y = {"station:o", title = nothing},

            color = {
                "data_fraction:q",
                aggregate = "sum",
                title = "Noise Data Fraction",
                scale = {domain= {unionWith = [0.0, 1.0]},
                         # scheme = "yellowgreenblue" # colormap of data fraction
                         scheme = "redyellowgreen" # colormap of data fraction
                        },
                legend={
                        gradientLength=300, # legend hight
                        gradientThickness=20, # legend width
                        labelFontSize=16}
            },

            config = {axis = {grid = true}},
        ))



    elseif  lowercase(xtimeunit)=="month"
        println("test")
        p = (df_all |> @vlplot(
            width = figsize[1],
            height = figsize[2],
            mark = :rect,
            x = {
                field = :date,
                title = "Time",
                timeUnit = :datemonthyearhoursminutesseconds,
                axis = {
                    format = xtimeaxisformat,
                    labelAngle = xlabelrotation,
                    labelOverlap = false,
                    #===here is the custumize of x axis span==#
                    labelColor={
                        condition={
                            test={
                                field="value",
                                timeUnit=:date,
                                equal={date=1}
                            },
                            value="black"
                        },
                        value=nothing
                    },
                    tickColor={
                        condition={
                            test={
                                field="value",
                                timeUnit=:date,
                                equal={date=1}
                            },
                            value="black"
                        },
                        value=nothing
                    },
                    #================#
                },
            },

            y = {"station:o", title = nothing},

            color = {
                "data_fraction:q",
                aggregate = "sum",
                title = "Noise Data Fraction",
                scale = {domain= {unionWith = [0.0, 1.0]},
                         # scheme = "yellowgreenblue" # colormap of data fraction
                         scheme = "redyellowgreen" # colormap of data fraction
                        },
                legend={
                        gradientLength=300, # legend hight
                        gradientThickness=20, # legend width
                        labelFontSize=16}
            },

            config = {axis = {grid = true}},
        ))

    else
        p = (df_all |> @vlplot(
            width = figsize[1],
            height = figsize[2],
            mark = :rect,
            x = {
                field = :date,
                title = "Time",
                timeUnit = :datemonthyearhoursminutesseconds,
                axis = {
                    format = xtimeaxisformat,
                    labelAngle = xlabelrotation,
                    labelOverlap = false,
                },
            },

            y = {"station:o", title = nothing},

            color = {
                "data_fraction:q",
                aggregate = "sum",
                title = "Noise Data Fraction",
                scale = {domain= {unionWith = [0.0, 1.0]},
                         # scheme = "yellowgreenblue" # colormap of data fraction
                         scheme = "redyellowgreen" # colormap of data fraction
                        },
                legend={
                        gradientLength=300, # legend hight
                        gradientThickness=20, # legend width
                        labelFontSize=16}
            },

            config = {axis = {grid = true}},
        ))
    end

    return p

end
