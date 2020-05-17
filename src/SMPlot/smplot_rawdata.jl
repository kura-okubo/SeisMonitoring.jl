using Plots
using SeisMonitoring: assemble_seisdata, scan_stations
"""
    smplot_rawdata(filename::String, station::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    dpi::Int=100)

Plot raw data within starttime and endtime at station. All channels are plotted.

# Argument

- `filename::String`                  : raw data (or RemovedEQ data) jld2 filename saved in SeisMonitoring.jl format.
- `fodir::String`                     : figure output directory.
- `station::String`                   : station name to plot.
- `starttime::Union{String, DateTime}`: starttime to plot.
- `endtime::Union{String, DateTime}`  : endtime plot.
- `max_point_num::Int=1e6`            : maximum number point for plotting in order to avoid too many sampling points during plot.
- `dpi::Int=100`                      : output figure dpi.
"""
function smplot_rawdata(filename::String, fodir::String, station::String,
                        starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                        max_point_num::Int=100000, dpi::Int=100) #NOTE:Int=1e6 is not working.

    ispath(filename) ? (fi = jldopen(filename, "r")) : error("$(filename) not found.")
    typeof(starttime) == String && (starttime = DateTime(starttime))
    typeof(endtime) == String && (endtime = DateTime(endtime))

    StationDict = scan_stations(filename)

    !haskey(StationDict, station) && error("$(station) is not saved in $(filename). Please check the JLD2 file. ")

    S = SeisData()

    for netstachan in StationDict[station]
            println("read $(netstachan)")
            Ch = assemble_seisdata(netstachan, fi, starttime, endtime)
            isnothing(Ch) || isempty(Ch) && continue
            push!(S, Ch)
    end

    isempty(S) && error("no data found at $(station) $(string(starttime)) - $(string(endtime)).")

    # plot traces
    plot(layout = (S.n, 1))
    layout = (S.n, 1)
    p_all = []
    for i =1:S.n
        S1 = S[i]
        # compute plot span
        Npts = length(S1.x)
        @show plotspan = ceil(Int, Npts/max_point_num)
        tvec_unix = S1.t[1,2]*1e-6 .+ (collect(1:plotspan:Npts).-1) ./ S1.fs # first sampling point is at S1.t[1,2]

        if i<S.n
            xformatter=x->""
        else
            xformatter=x->Dates.format.(Dates.unix2datetime.(x),"yyyy/m/d HH:MM")
        end

        Plots.plot!(tvec_unix, S1.x[1:plotspan:Npts], xformatter=xformatter,
               linewidth = 1.0, linecolor=:black, title=S1.id, xrotation = 0,
               xtickfontsize=8,ytickfontsize=8, xguidefontsize=8, yguidefontsize=8,
               titlefontsize=8, label="", subplot=i)
        Plots.ylabel!(S1.units)
    end

    p = Plots.plot!(link=:x, dpi=dpi)
    figname = joinpath(fodir, join([station, string(starttime),  string(endtime)], "__"))
    Plots.savefig(p, figname)
    close(fi)
end
