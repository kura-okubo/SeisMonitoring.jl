using SeisIO, SeisNoise, Dates, DSP, Statistics, StatsBase, LinearAlgebra, ColorSchemes, JLD2, Plots
using SeisMonitoring: assemble_seisdata, smooth_withfiltfilt
"""
    compute_psdpdf(S::SeisChannel, starttime::DateTime, endtime::DateTime)

Return mean and probablistic density function of seismic noise periodgram.
"""
function compute_psdpdf(S::SeisChannel, figname::String, segments_length::Real;
        figsize=(1200, 800), fmt::String="png", xlims=(0.01, 5.0), ylims=(-200, -80), clims=(0.0, 0.3))

    println("-------Plot Periodogram--------")
    # 1. truncate with respect to the segments
    R = RawData(S, segments_length, segments_length)
    # 2. Preprocessing
    detrend!(R)
    taper!(R)
    # 3. Compute Periodgram
    P = DSP.periodogram(R.x, fs=R.fs)
    # 4. Apply  smoothing NOTE:here using filtfilt instead of SeisNoise.smooth()
    P.power[:, :] = smooth_withfiltfilt(P.power, window_len=7, window=:rect) #immutable struct of type Periodogram2 cannot be changed
    # 5. Plot periodogram
    p1=plot(size=figsize)

    tracespan = ceil(Int, size(P.power, 2)/1000) # maximum 1000 traces

    for i = 1:tracespan:size(P.power, 2)
        p1 = plot!(P.freq1, DSP.pow2db.(P.power[:,i]),
            xlims=xlims, ylims=ylims, xscale=:log10, label="", color=:gray, alpha=0.4)
    end
    # plot median
    p1 = plot!(P.freq1, DSP.pow2db.(mean(P.power, dims=2)), label="median", color=:red)
    # min and max
    p1 = plot!(P.freq1, DSP.pow2db.(minimum(P.power, dims=2)), label="minimum", color=:blue)
    p1 = plot!(P.freq1, DSP.pow2db.(maximum(P.power, dims=2)), label="maximum", color=:purple)
    # plot waterlevel guide
    freqrange = [minimum(xlims); maximum(xlims)]
    for water_level in [1e-3, 1e-4, 1e-5]
        # compute waterlevel for deconvolution
        wl = DSP.pow2db(water_level * mean(mean(P.power, dims=2),dims=1)[1])
        p1 = plot!(freqrange, [wl; wl], label="water_level = $(water_level)", color=:black, linestyle=:dot)
    end

    s_str = string(u2d(S.t[1,2]*1e-6))[1:19]
    et = S.t[1,2]*1e-6 + (S.t[end,1]-1)/S.fs
    e_str=string(u2d(et))[1:19]
    title!("$(S.id) $(s_str) - $(e_str)")
    xlabel!("Frequency [Hz]")
    ylabel!("Power[10log10 (m**2/sec**2/Hz)][dB]")
    savefig(p1, "periodogram_$(figname).$(fmt)")

    # compute probability density function of power spectral density (PDFPSD)
    # following McNamara and Buland (2004), using 1-db binned  between -200:-80 db histogram to obtain pdf.
    println("-------Plot PDFPSD--------")

    bins = -200:1.0:-80
    rfft_halflen =  round(Int, size(P.power,1)/2) # assuming real signal, compute half of the fft.
    ppsd = Array{Float64, 2}(undef, length(bins)-1, rfft_halflen) #(number of bins, frequency)

    for i =1:rfft_halflen
         h = StatsBase.fit(Histogram, DSP.pow2db.(P.power[i,:]), bins)
         hn = normalize(h, mode=:pdf)
         ppsd[:, i] = hn.weights
    end

    # compute logspace bins in frequency to better plot contourf of probability density function
    Nplotfreq = 1000
    xq = 10 .^(range(log10(xlims[1]), stop=log10(xlims[2]), length=Nplotfreq))
    ixq = []
    for i = 1:length(xq)
        push!(ixq, findfirst(x -> x >= xq[i], P.freq1))
    end
    unique!(ixq)

    # load pqlx colormap
    c_pqlx = jldopen("pqlx.jld2", "r") do fi; fi["pqlx.colors"]; end
    loadcolorscheme(:pqlx, c_pqlx, "pqlx color", "for pdfpsd plot")
    #===#

    # NOTE: heatmap is not clean to plot
    # heatmap(P.freq1[ixq], collect(bins)[1:end-1], ppsd[:, ixq],
    #         xscale=:log10, clims=clims, xlims=(0.01,5.0), ylims=(-200,-80),
    #         color=heatmapcolor)

    contourf(P.freq1[ixq], collect(bins)[1:end-1], ppsd[:, ixq], color=:pqlx,
    levels=30, clims=clims, xscale=:log10, xlims=xlims, ylims=ylims, linewidth=0.0,
    size=figsize, colorbar_title="Probability")

    title!("$(S.id) $(s_str) - $(e_str)")
    xlabel!("Frequency [Hz]")
    ylabel!("Power[10log10 (m**2/sec**2/Hz)][dB]")
    savefig("pdfpsd_$(figname).$(fmt)")
end


"""
    smplot_pdfpsd(filename::String, figname::String, station::String, starttime::DateTime,
        endtime::DateTime; segments_length::Real=3600, figsize=(1200, 800), xlims = (0.01, 5.0),
        ylims = (-200, -80), clims = (0.0, 0.3), fmt = "png")

Plot the probability density function of power spectral density (PSDPDF) following McNamara and Buland (2004).

# Argument
-`filename::String`: absolute/relative path to SeisMonitoring.jl format data. (e.g. "OUTPUT/RawData.jld2")
-`figname::String`: figure names. (e.g. "fig1" produces "periodgram_fig1.fmt" and "pdfpsd_fig1.fmt")
-`channel::String`: station, channel name (e.g. BP.CCRB..BP1)
-`starttime::DateTime`: start time of continuous data.
-`endtime::DateTime`: end time of continuous data.
-`segments_length::Real=3600`: [second] segments length to compute power spectral density. 3600sec is used in McNamara and Buland (2004).

# Option
-`figsize=(1200, 800)`: figure size to be plotted.
-`xlims = (0.01, 5.0)`: xlims
-`ylims = (-200, -80)`: ylims
-`clims = (0.0, 0.3)`: color range of pdfpsd.
-`fmt = "png"`: figure format
"""
function smplot_pdfpsd(filename::String, figname::String, channel::String, starttime::DateTime,
    endtime::DateTime, segments_length::Real=3600; figsize=(1200, 800), xlims = (0.01, 5.0),
    ylims = (-200, -80), clims = (0.0, 0.3), fmt = "png")

    S1 = assemble_seisdata(channel, fi, starttime, endtime, data_contents_fraction=0.0)

    compute_psdpdf(S1, figname, segments_length, figsize=figsize, fmt=fmt, xlims=xlims, ylims=ylims, clims=clims)
end

# # S = rseis("psdtest.seisio")[1]
# filename = "/Volumes/Kurama_20190821/kurama/research/SeisMonitoring_dev/testproject_2_OUTPUT/seismicdata/EQRemovedData.jld2"
# # channel = "BP.EADB..BP1"
# fi = jldopen(filename, "r");
# starttime = DateTime("2004-04-01")
# endtime = DateTime("2004-04-07")
# segments_length = 3600 #time series length (e.g. 1 hour for McNamara and Buland, 2004)
#
# # figname = channel
# fmt = "png"
#
# figsize=(800, 600)
# xlims = (0.01, 5.0)
# ylims = (-200, -80)
# clims = (0.0, 0.3)
#
# smplot_pdfpsd(filename, "BP.CCRB..BP1", "BP.CCRB..BP1", starttime, endtime, segments_length,
#             figsize=figsize, xlims=xlims,ylims=ylims,clims=clims,fmt=fmt)
#
# smplot_pdfpsd(filename, "BP.EADB..BP1", "BP.EADB..BP1", starttime, endtime, segments_length,
#             figsize=figsize, xlims=xlims,ylims=ylims,clims=clims,fmt=fmt)
#
# smplot_pdfpsd(filename, "BP.FROB..BP1", "BP.FROB..BP1", starttime, endtime, segments_length,
#             figsize=figsize, xlims=xlims,ylims=ylims,clims=clims,fmt=fmt)
#



#===NOTE: if the pdlx.jld2 is missing, run the script below to make colorschemes from (r,g,b) text file===#
# lines = readlines("pqlx.txt", keep=true)
# c_pqlx = []
# for line in lines
#     r, g, b = parse.(Float64, split(split(line)[1], ","))
#     push!(c_pqlx, (r, g, b))
# end
# pqlx = ColorScheme([RGB{Float64}(c_pqlx[i][1],  c_pqlx[i][2], c_pqlx[i][3] ) for i in 1:length(c_pqlx)])
# @save "pqlx.jld2" pqlx.colors
# #===========================================================================#
