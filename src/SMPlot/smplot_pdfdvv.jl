using SeisIO, SeisNoise, JLD2, Dates, DataFrames, CSV, StatsBase, Statistics, ColorSchemes, LinearAlgebra, Plots

"""
    smplot_pdfdvv()

Plot the statistics of change of dvv history.

---

# NOTE

Please run first `smstats_read()` to make `monitoring_stats.csv`. `smplot_pdfdvv()` assemble
dvv history database.

---


"""
function smplot_pdfdvv(statsfile::String, fodir::String, starttime::DateTime, endtime::DateTime, time_bin_length::Real;
    cc_threshold::Float64=0.6, network_option=["all"], compontents_option::AbstractArray=["all"],
    plot_maxdvv::Float64=0.05, number_of_dvvbins::Int = 30, minimum_paircount::Int=1, plottimeunit::String="day",
    figsize = (1200, 600), clims=(0.0, 0.3), xlims = [], ylims = (-0.03, 0.03), xrotation::Real=-45,
    fmt="png")

    #1. read dvvfile
    df = CSV.read(statsfile)

    stbins = range(d2u(starttime), stop = (d2u(endtime)-time_bin_length), step = time_bin_length)
    etbins = stbins .+ time_bin_length
    mtbins = (stbins .+ etbins)/2

    # scan freqband
    freqbands = sort(unique(df.freqband))

    # loop for each freqband

    DvvDicts = []

    for freqband in freqbands

        # filter with freqband
        df_filtered = filter(:freqband => x -> (x == freqband) , df)

        # filtering with network option,
        df_filtered_net = similar(df_filtered, 0)
        for net in network_option
            println("$(net)-$(net)")
            # if network is chosen, select ccs only within the network, while "all" computes across networks
            df_temp = filter(:networks => x -> (x == "$(net)-$(net)" || "all" ∈ network_option) , df_filtered)
            append!(df_filtered_net, df_temp)
        end
        #rename for later proces
        df_filtered = df_filtered_net

        # filtering with component option
        df_filtered = filter(:components => x -> (x ∈ compontents_option || "all" ∈ compontents_option) , df_filtered)

        # filtering with cc_threshold
        df_filtered = filter(:cc_dvv => x -> (x >= cc_threshold) , df_filtered)

        dvvbins = range(-plot_maxdvv, stop=plot_maxdvv, length=number_of_dvvbins+1)

        # compute count, psd and mean with respect to time bins
        DvvDict = Dict(
            "freqband" => freqband,
            "T" => DateTime[], #DateTime
            "dvv_mean" => Union{Float64, Missing}[],
            "pdfdvv" => Array{Union{Float64, Missing}, 2}(undef, number_of_dvvbins, 0),
            "count_pairs" => Int[],
        )


        for (i, mtbin) in enumerate(mtbins)
            stbin = u2d(stbins[i])
            etbin = u2d(etbins[i])

            println((stbin, etbin))

            push!(DvvDict["T"], u2d(mtbin))

            global df_binned = filter(:date => x -> (stbin <= x < etbin) , df_filtered)
            println(df_binned)
            dvv_all = df_binned.dvv

            paircount = length(dvv_all)

            if paircount < minimum_paircount
                # this time bin does not have enough station pairs above threshold
                # append to DvvDict
                push!(DvvDict["dvv_mean"], missing)
                DvvDict["pdfdvv"] = hcat(DvvDict["pdfdvv"], fill(missing, number_of_dvvbins))
                push!(DvvDict["count_pairs"], paircount)
                continue;
            end

            #compute pdf of dvv
             h = StatsBase.fit(Histogram, dvv_all, dvvbins)
             hn = normalize(h, mode=:probability)
             println(hn.weights)
             # append to DvvDict
             push!(DvvDict["dvv_mean"], Statistics.mean(dvv_all))
             DvvDict["pdfdvv"] = hcat(DvvDict["pdfdvv"], hn.weights)
             push!(DvvDict["count_pairs"], paircount)
        end

        # plot pdf-mean of dvv and bar count.
        plot(layout = grid(2, 2, heights=[0.6, 0.4], widths=[0.9, 0.1]), link=:x)

        # plot title
        figtitle = "$(string(starttime)) - $(string(endtime))  $(DvvDict["freqband"])Hz"
        # plot!(annotation=(0.5,0.5, figtitle), framestyle = :none, subplot=1, fontsize=12)

        # 1. pdf and mean
        if lowercase(plottimeunit) == "year"
            xaxisformat = "yyyy"
            xtickid = findall(x -> (Dates.Month(u2d(x)).value == 1 && Dates.Day(u2d(x)).value == 1), mtbins)
            !isempty(xtickid) ? xticks = mtbins[xtickid] :  (@warn("no xticks with $(plottimeunit) plottimeunit. plot auto."); xticks=:auto)
        elseif lowercase(plottimeunit) == "month"
            xaxisformat = "yyyy-m"
            xtickid = findall(x -> (Dates.Day(u2d(x)).value == 1), mtbins)
            !isempty(xtickid) ? xticks = mtbins[xtickid] :  (@warn("no xticks with $(plottimeunit) plottimeunit. plot auto."); xticks=:auto)
        elseif lowercase(plottimeunit) == "day"
            xaxisformat = "yyyy-m-d"
            xticks=:auto
        end


        xformatter=x->Dates.format.(Dates.unix2datetime.(x), xaxisformat)

        isempty(xlims) &&  (xlims = (d2u(starttime), d2u(endtime)))

        # load colormap
        smplot = jldopen(joinpath(splitdir(pathof(SeisMonitoring))[1], "SMPlot/smplot_pararainbow.jld2"), "r") do fi; fi["smplot.colors"]; end
        loadcolorscheme(:smplot_pararainbow, smplot, "smplot color", "for smplot")

        plot!(link=:x, subplot=1)
        contourf!(mtbins, collect(dvvbins)[1:end-1], DvvDict["pdfdvv"], color=:smplot_pararainbow,
        levels=30, clims=clims, xformatter=xformatter, xlims=xlims, ylims=ylims, linewidth=0.0,
        size=figsize, colorbar_title="Probability", subplot=1, label="", ylabel="dv/v",
        xrotation=xrotation, frame=:box, xticks = xticks,
        title = figtitle, colorbar=false)

        plot!(mtbins, DvvDict["dvv_mean"], size=figsize, color=:magenta, subplot=1, linewidth = 3.0,
         xformatter=xformatter, xticks = xticks, label = "mean", link=:x)

        # plot only colorbar to share the x axis between contour and bar.
        # NOTE: Currently there is no way to plot only colorbar with julia.
        heatmap!(fill(missing,1,1), color=:smplot_pararainbow, levels=30, clims=clims, colorbar_title="Probability",
         subplot=2, label="",colorbar=true, framestyle = :none)
        plot!(grid = false, showaxis = false, subplot=4)

        # 2. bar of pairs count
        bar!(mtbins, DvvDict["count_pairs"], fillcolor=:blue, fillalpha=0.8,
            xrotation=xrotation, frame=:box, xticks = xticks,
            xformatter=xformatter, label="", ylabel="count", subplot=3, link=:x)

        xlims!((d2u(starttime), d2u(endtime))) # NOTE: This is the way to share the axis
        p = Plots.plot!(link=:x)

        # save figure
        figname = joinpath(fodir, join(["psfdvv_", string(starttime), string(endtime), freqband*"Hz"], "__"))*".$(fmt)"
        Plots.savefig(p, figname)
        # save JLD2 file to reproduce the figure
        jldopen("$(figname).jld2", "w") do fo
            fo["DvvDict"] = DvvDict
        end
        push!(DvvDicts, DvvDict)
    end

    return DvvDicts
end
