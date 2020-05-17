using Plots

"""
    smplot_shorttimestack(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                    channelpair::AbstractArray=["all"], dpi::Int=100)

Plot Stacked corrdata within the requested start and endtime.
"""
function smplot_shorttimestack(filename::String, fodir::String, starttime::Union{String, DateTime}, endtime::Union{String, DateTime};
                channelpair::AbstractArray=["all"], dpi::Int=100)

    ispath(filename) ? fi = jldopen(filename, "r") : error("$(filename) not found.")
    typeof(starttime) == String && (starttime = DateTime(starttime))
    typeof(endtime) == String && (endtime = DateTime(endtime))

    for stachankey = keys(fi)
        sta1, sta2 = split(stachankey, "-")
        @show comp = sta1[end]*sta2[end]
        if comp ∈ channelpair || "all" ∈ channelpair
            stachan = fi[stachankey]
            for tkey in keys(stachan)
                st, et = DateTime.(split(tkey, "--"))
                if starttime <= st && et <= endtime
                # plot all freq within this timewindow

                   p = plot(bg=:white)
                   stack_method = ""

                   # sort by frequency
                   for (ifreq, freqkey) = enumerate(sort(keys(stachan[tkey])))
                       C1 = stachan[joinpath(tkey, freqkey)]
                       # plot coda_window box
                       fillbox = C1.misc["fillbox"]
                       yshift = 2*float(ifreq-1)
                       yshift_box = [yshift, yshift]
                       plot!(fillbox[1:2], yshift_box,
                           fillrange=[yshift_box.-0.9], fillalpha=0.1, c=:orange,
                           label="", linealpha=0.0)
                       plot!(fillbox[1:2], yshift_box,
                           fillrange=[yshift_box.+0.9], fillalpha=0.1, c=:orange,
                           label="", linealpha=0.0)
                       plot!(fillbox[3:4], yshift_box,
                           fillrange=[yshift_box.-0.9], fillalpha=0.1, c=:orange,
                           label="", linealpha=0.0)
                       plot!(fillbox[3:4], yshift_box,
                           fillrange=[yshift_box.+0.9], fillalpha=0.1, c=:orange,
                           label="", linealpha=0.0)

                       if haskey(C1.misc, "reference")
                           normamp = maximum(abs.(C1.misc["reference"])) # normalized amplitude by reference
                           Plots.plot!(C1.misc["timelag"], C1.misc["reference"]./normamp .+ yshift, color=:black, label="")
                       else
                           normamp = maximum(C.corr) # normalized amplitude by current corr
                       end

                       Plots.plot!(C1.misc["timelag"], C1.corr./normamp .+ yshift, color=ifreq, label=freqkey*"[Hz]")
                       if ifreq == 1
                           max_coda_length = maximum(abs.(fillbox))
                           Plots.xlims!((-1.2*max_coda_length, 1.2*max_coda_length))
                       end

                       isempty(stack_method) && (stack_method = "$(C1.misc["stack_method"])-$(C1.corr_type)")
                   end



                   # plot stretching time window
                   fontsize = 8
                   p = Plots.plot!(xlabel = "Time lag [s]",  ylabel = "NCFs", yformatter=_->"",
                       xtickfontsize=fontsize,ytickfontsize=fontsize,
                       xguidefontsize=fontsize,yguidefontsize=fontsize,legendfontsize=fontsize,
                       legend=:outertopright, title = "$(stachankey): $(tkey): stack method:$(stack_method)\nBlack: reference, Color: shorttime stack",
                       titlefontsize=8, dpi=dpi)

                   figname = joinpath(fodir, join(["stack", stachankey, tkey], "_"))
                   Plots.savefig(p, figname)
                end
            end
        end
    end
    close(fi)
end
