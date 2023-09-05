using SeisMonitoring: assemble_corrdata, cc_medianmute!
using SeisIO, SeisNoise, JLD2, Plots, Dates
"""
    smplot_corrdata(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    channelpair::AbstractArray=["all"], dpi::Int=100)

Plot CorrData within the requested start and endtime.

# NOTE
    The `linear` stacked trace is shown at bottom panel. Please use smplot_stackcc() to plot the stacked data by SeisMonitoring.seisstack().
"""
function smplot_corrdata(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                channelpair::AbstractArray=["all"], dpi::Int=100, plotmaxlag::Float64=0.0,
                cc_medianmute_max::Float64=2.0, cc_medianmute_min::Float64=0.0, MAX_MEM_USE::AbstractFloat=4.0)

    ispath(filename) ? fi = jldopen(filename, "r") : error("$(filename) not found.")
    typeof(starttime) == String && (starttime = DateTime(starttime))
    typeof(endtime) == String && (endtime = DateTime(endtime))

    # for stachankey = keys(fi)
    # println(stachankey)
    stachanpair = splitdir(filename)[2][1:end-5] # BP.LCCB-BP.VCAB-11
    sta1, sta2, comp = split(stachanpair, "-")

    # sta1, sta2 = split(stachankey, "-")
    # comp = sta1[end]*sta2[end]

    if comp ∈ channelpair || "all" ∈ channelpair

        # get frequency band
        freqency_band = Float64[]
        # stachan = fi[stachankey]
        # for tkey in keys(stachan)
        for tkey in keys(fi)
            st, et = DateTime.(split(tkey, "--"))
            if starttime <= st && et <= endtime
            # plot all freq within this timewindow
                for freqkey in keys(fi[tkey])
                    freqmin, freqmax = parse.(Float64, split(freqkey, "-"))
                    freqmin ∉ freqency_band && push!(freqency_band, freqmin)
                    freqmax ∉ freqency_band && push!(freqency_band, freqmax)
                end
            end
        end

        isempty(freqency_band) && error("no corrdata within $(starttime) and $(endtime) is found.")
        sort!(freqency_band)

        # plot corrdata with respect to frequency band
        # println(C_all)

        # stack with respect to frequency band
        Nfreqband = length(freqency_band) - 1
        freqband = map(i -> [freqency_band[i], freqency_band[i+1]], 1:Nfreqband)

        for fb in freqband
            freqmin, freqmax = fb
            freqkey = join([string(freqmin), string(freqmax)], "-")

            # assemble corrdata on stachankey
            C1, CorrData_Buffer = assemble_corrdata(fi,starttime,endtime,freqkey,
                                                MAX_MEM_USE=MAX_MEM_USE, min_cc_datafraction=0.0)

            # C1 = C_all[freqkey]
            cc_medianmute!(C1, cc_medianmute_max, cc_medianmute_min)

            # modify plot maxlag
            if isempty(C1.t)
                @warn("empty corrdata. skip"); continue
            end
            tkey = join([starttime, endtime], "--")
            SeisNoise.corrplot(C1)
            p = Plots.plot!(title = "$(stachanpair): $(tkey): $(freqkey) Hz", titlefontsize=8, dpi=dpi)
            !iszero(plotmaxlag) && xlims!(-plotmaxlag, plotmaxlag)
            @show figname = joinpath(fodir, join([stachanpair, tkey, freqkey], "__"))*".png"
            Plots.savefig(p, figname)
        end
    end
    # end

    close(fi)

end



# fi = jldopen("/Volumes/Kurama_20200621/research/SeisMonitoring_BP/Parallelization_dev/init/BP_cc_parallel_OUTPUT/cc/BP.LCCB.40.SP1-BP.LCCB.40.SP1.jld2", "r")
# stachanpair = "BP.CCRB..BP1-BP.EADB..BP1"
# starttime = DateTime("2015-01-01T00:00:00")
# endtime = DateTime("2015-01-31T00:00:00")
# fodir = "./"
# filename = "/Volumes/Kurama_20200621/research/SeisMonitoring_BP/dev_assemblecorrdata_collected/BP.LCCB-BP.LCCB-13.jld2"
# smplot_corrdata(filename, fodir, starttime, endtime)
