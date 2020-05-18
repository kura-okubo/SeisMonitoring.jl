using SeisMonitoring: assemble_corrdata
"""
    smplot_corrdata(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    channelpair::AbstractArray=["all"], dpi::Int=100)

Plot CorrData within the requested start and endtime.

# NOTE
    The `linear` stacked trace is shown at bottom panel. Please use smplot_stackcc() to plot the stacked data by SeisMonitoring.seisstack().
"""
function smplot_corrdata(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                channelpair::AbstractArray=["all"], dpi::Int=100, MAX_MEM_USE::AbstractFloat=4.0)

    ispath(filename) ? fi = jldopen(filename, "r") : error("$(filename) not found.")
    typeof(starttime) == String && (starttime = DateTime(starttime))
    typeof(endtime) == String && (endtime = DateTime(endtime))

    for stachankey = keys(fi)
        println(stachankey)
        sta1, sta2 = split(stachankey, "-")
        comp = sta1[end]*sta2[end]
        if comp ∈ channelpair || "all" ∈ channelpair

            # get frequency band
            freqency_band = Float64[]
            stachan = fi[stachankey]
            for tkey in keys(stachan)
                st, et = DateTime.(split(tkey, "--"))
                if starttime <= st && et <= endtime
                # plot all freq within this timewindow
                    for freqkey in keys(stachan[tkey])
                        freqmin, freqmax = parse.(Float64, split(freqkey, "-"))
                        freqmin ∉ freqency_band && push!(freqency_band, freqmin)
                        freqmax ∉ freqency_band && push!(freqency_band, freqmax)
                    end
                end
            end

            isempty(freqency_band) && error("no corrdata within $(starttime) and $(endtime) is found.")
            sort!(freqency_band)

            # assemble corrdata on stachankey
            C_all, CorrData_Buffer = assemble_corrdata(fi,stachankey,starttime,endtime,freqency_band,
                                                MAX_MEM_USE=MAX_MEM_USE, min_cc_datafraction=0.0)

            # plot corrdata with respect to frequency band
            println(C_all)
            for freqkey in keys(C_all)
                C1 = C_all[freqkey]
                if isempty(C1.t)
                    @warn("empty corrdata. skip"); continue
                end
                tkey = join([starttime, endtime], "--")
                SeisNoise.corrplot(C1)
                p = Plots.plot!(title = "$(stachankey): $(tkey): $(freqkey) Hz", titlefontsize=8, dpi=dpi)
                @show figname = joinpath(fodir, join([stachankey, tkey, freqkey], "__"))*".png"
                Plots.savefig(p, figname)
            end
        end
    end

    close(fi)

end
