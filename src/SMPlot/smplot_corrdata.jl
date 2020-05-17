using Plots
"""
    smplot_corrdata(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    channelpair::AbstractArray=["all"], dpi::Int=100)

Plot CorrData within the requested start and endtime.

# NOTE
    The `linear` stacked trace is shown at bottom panel. Please use smplot_stackcc() to plot the stacked data by SeisMonitoring.seisstack().
"""
function smplot_corrdata(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                channelpair::AbstractArray=["all"], dpi::Int=100)

    ispath(filename) ? fi = jldopen(filename, "r") : error("$(filename) not found.")
    typeof(starttime) == String && (starttime = DateTime(starttime))
    typeof(endtime) == String && (endtime = DateTime(endtime))

    for stachankey = keys(fi)
        println(stachankey)
        sta1, sta2 = split(stachankey, "-")
        comp = sta1[end]*sta2[end]
        if comp ∈ channelpair || "all" ∈ channelpair
            stachan = fi[stachankey]
            for tkey in keys(stachan)
                st, et = DateTime.(split(tkey, "--"))
                if starttime <= st && et <= endtime
                # plot all freq within this timewindow
                    for freqkey in keys(stachan[tkey])
                        C1 = stachan[joinpath(tkey, freqkey)]
                        SeisNoise.corrplot(C1)
                        p = Plots.plot!(title = "$(stachankey): $(tkey): $(freqkey) Hz", titlefontsize=8, dpi=dpi)
                        figname = joinpath(fodir, join([stachankey, tkey, freqkey], "__"))
                        Plots.savefig(p, figname)
                    end
                end
            end
        end
    end

    close(fi)

end
