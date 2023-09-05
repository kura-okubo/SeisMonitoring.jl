using SeisIO, SeisNoise, Dates, DSP, Statistics, StatsBase, LinearAlgebra, ColorSchemes, JLD2, Plots, DataFrames, CSV, Measures
using SeisMonitoring: assemble_seisdata, smooth_withfiltfilt

"""
    compute_psdpdf(S::SeisChannel, fodir::String, segments_length::Real;
	        overlap::Real=segments_length/2, smoothing_window_len::Real=7, figsize=(1200, 800),
			fmt::String="png", xlims=(0.01, 5.0), ylims=(-200, -80), clims=(0.0, 0.3))

Compute power spectrum density - probablistic density function.
Refering process flow from obspy.signal.spectral_estimation.PPSD: https://docs.obspy.org/_modules/obspy/signal/spectral_estimation.html#PPSD

# Argument
-`S::SeisChannel`		: SeisChannel containing the waveform to compute pdfpsd.
-`fodir::String`		: figure output directory
-`segments_length::Real`: [s] Length of segment to compute windowed psd.

# Option
-`overlap::Real`		: [s] Length of overlap between segments to compute windowed psd.
-`smoothing_window_len::Real` : [npts] Smoothing points.
-`figsize=(1200, 800)`: figure size to be plotted.
-`xlims = (0.01, 5.0)`: xlims
-`ylims = (-200, -80)`: ylims
-`clims = (0.0, 0.3)`: color range of pdfpsd.
-`fmt = "png"`: figure format

# Output
-`freqvec`	: [Hz] Vector of frequency for the ppsd.
-`bins`		: [dB] bins of amplitude for the ppsd.
-`ppsd`		: [probability] probability density function of psd.
-`(percentile_10, percentile_50, percentile_90)`		: Vector of 10th, 50th, 90th percentiles.

# Note
About instrumental response removal, the obspy.signal.spectral_estimation.PPSD do the removal during the computation of PPSD.
i.e. the raw data (tr.data)


Updated on 2021.05.25 Kurama Okubo
"""
function compute_psdpdf(S::SeisChannel, fodir::String, segments_length::Real;
        overlap::Real=segments_length/2, smoothing_window_len::Real=7, figsize=(1200, 800),
		fmt::String="png", xlims=(0.01, 5.0), ylims=(-200, -80), clims=(0.0, 0.3), eps::Real=1e-24)

	println("-------Compute PPSD $(S.id)--------")

	!ispath(fodir) && mkpath(fodir)

	R1 = RawData(S, segments_length, overlap)
	Nseg = size(R1.x, 2)
	Nwin = trunc(Int, segments_length*R1.fs / 4) # to make 13 segments on each trace with 75 % overlap, following obspy and McNamara and Buland (2004)
	Noverlap = trunc(Int, 0.75 * Nwin)
	tukey10(x) = tukey(x, 0.2) # see obspy fft_taper, which uses 10% cosine taper following McNamara and Buland (2004)

	# compute dummy psd for freq vector
	Pseg_dum = welch_pgram(R1.x[:, 1], Nwin, Noverlap, onesided=true, fs=R1.fs, window=tukey10)
	freqvec = Pseg_dum.freq
	ω = 2.0 * π .* freqvec
	bins = -200:1.0:-80 # Fix bins for simplicity

	p_all = zeros(Float64, length(Pseg_dum.power), Nseg)
	ppsd =  zeros(Float64, length(bins)-1, length(freqvec)) #(number of bins, frequency)

	# compute periodogram using welch method
	t2 = @elapsed for i in 1:Nseg
		Pseg = welch_pgram(R1.x[:, i], Nwin, Noverlap, onesided=true, fs=R1.fs, window=tukey10)
		#NOTE: Convert from velocity to acceleration power
		if S.units == "m/s"
			p_all[:, i] = Pseg.power .* ω.^2 # note that this is power/Hz normalized by window when using oneside=true
		elseif S.units == "m/s/s"
			p_all[:, i] = Pseg.power
		elseif S.units == "m"
			p_all[:, i] = Pseg.power .* ω.^4 # displacement -> acceleration
		else
			error("Unit of data $(S.units) is unknown.")
		end
	end

	p_all[:, :] = smooth_withfiltfilt(p_all, window_len=smoothing_window_len, window=:rect) #immutable struct of type Periodogram2 cannot be changed

	# Compute 10th, 50th and 90th percentiles
	percentile_10 = zeros(length(freqvec), 1)
	percentile_50 = zeros(length(freqvec), 1)
	percentile_90 = zeros(length(freqvec), 1)

	for i =1:length(freqvec)
		p10 = StatsBase.percentile(p_all[i, :], 10)
		p50 = StatsBase.percentile(p_all[i, :], 50)
		p90 = StatsBase.percentile(p_all[i, :], 90)
		percentile_10[i] = p10 > eps ? DSP.pow2db(p10) : NaN
		percentile_50[i] = p50 > eps ? DSP.pow2db(p50) : NaN
		percentile_90[i] = p90 > eps ? DSP.pow2db(p90) : NaN
		# check if 50th percentile == median
		if ~isnan(percentile_50[i])
			@assert percentile_50[i] == DSP.pow2db.(Statistics.median(p_all[i, :]))
		end
	end

	# compute histogram on each frequency to obtain pdf
	for i =1:length(freqvec)
		psd_tmp = p_all[i,:]
		filter!(x -> x > eps, psd_tmp) # filter out the too small powers
		h = StatsBase.fit(Histogram, DSP.pow2db.(psd_tmp), bins)
		hn = normalize(h, mode=:pdf)
		ppsd[:, i] = hn.weights
	end

	# plot(bins[1:end-1], ppsd[:, 1:100:end])

	# test plot
	xlims=(0.01, 5.0)
	ylims=(-200, -50)

	Nplotfreq = 100
	xq = 10 .^(range(log10(xlims[1]), stop=log10(xlims[2]), length=Nplotfreq)) # logalithmically uniform octave space
	ixq = []
	for i = 1:length(xq)
	push!(ixq, findfirst(x -> x >= xq[i], freqvec))
	end
	unique!(ixq)

	# load pqlx colormap
	c_pqlx = jldopen(joinpath(splitdir(pathof(SeisMonitoring))[1], "SMPlot/smplot_pqlx.jld2"), "r") do fi; fi["pqlx.colors"]; end
	~haskey(colorschemes, :pqlx) && loadcolorscheme(:pqlx, c_pqlx, "pqlx color", "for pdfpsd plot");
	#===#

	t3 = @elapsed heatmap(freqvec[ixq], bins[1:end-1], ppsd[:, ixq], clims=(0.0, 0.3), color=:pqlx,
	 xlims=xlims, ylims=ylims, framestyle = :box, margin = 15mm, grid=true,
	 xticks=([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3], [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]),
	levels=30, xscale=:log10, linewidth=0.0, size=(800,600), colorbar_title="Probability",
	)

	#plot percentile
	plot!(freqvec[ixq], percentile_10[ixq], linestyle=:solid, linewidth = 2,linecolor = :cyan, label="10th%", legend=:topright)
	plot!(freqvec[ixq], percentile_50[ixq], linestyle=:solid, linewidth = 2,linecolor = :red, label="median")
	plot!(freqvec[ixq], percentile_90[ixq], linestyle=:solid, linewidth = 2,linecolor = :orange, label="90th%")

	s_str = string(u2d(S.t[1,2]*1e-6))[1:19]
	et = S.t[1,2]*1e-6 + (S.t[end,1]-1)/S.fs
	e_str=string(u2d(et))[1:19]
	title!("$(S.id) $(s_str) - $(e_str)")
	xlabel!("Frequency [Hz]")
	ylabel!("Power[10log10 (m**2/sec**4/Hz)][dB]")
	savefig(fodir*"/pdfpsd__$(S.id)__$(s_str)__$(e_str).$(fmt)")
	println(t2, t3)
    return freqvec, bins, ppsd, ixq, (percentile_10, percentile_50, percentile_90)

end


"""
    smplot_pdfpsd(fidir::String, fodir::String, station::String, starttime::DateTime,
        endtime::DateTime; segments_length::Real=3600, figsize=(1200, 800), xlims = (0.01, 5.0),
        ylims = (-200, -80), clims = (0.0, 0.3), fmt = "png")

Plot the probability density function of power spectral density (PSDPDF) following McNamara and Buland (2004).

# Argument
-`fidir::String`: absolute/relative path to SeisMonitoring.jl format data. (e.g. "OUTPUT/RawData.jld2")
-`fodir::String`:
-`netstalocchan::String`: station, channel name (e.g. BP.CCRB..BP1)
-`starttime::DateTime`: start time of continuous data.
-`endtime::DateTime`: end time of continuous data.
-`segments_length::Real=3600`: [second] segments length to compute power spectral density. 3600sec is used in McNamara and Buland (2004).


# Argument
-`S::SeisChannel`		: SeisChannel containing the waveform to compute pdfpsd.
-`fodir::String`		: figure output directory
-`segments_length::Real`: [s] Length of segment to compute windowed psd.

# Option
-`overlap::Real`		: [s] Length of overlap between segments to compute windowed psd.
-`smoothing_window_len::Real` : [npts] Smoothing points.
-`figsize=(1200, 800)`: figure size to be plotted.
-`xlims = (0.01, 5.0)`: xlims
-`ylims = (-200, -80)`: ylims
-`clims = (0.0, 0.3)`: color range of pdfpsd.
-`fmt = "png"`: figure format
"""
function smplot_pdfpsd(fidir::String, fodir::String, netstalocchan::String, starttime::DateTime,
    endtime::DateTime, segments_length::Real=3600; overlap::Real=segments_length/2, smoothing_window_len::Real=7, figsize=(1200, 800), xlims = (0.01, 5.0),
    ylims = (-200, -80), clims = (0.0, 0.3), fmt = "png")

    # fi = jldopen(fidir, "r")
    t1 = @elapsed S1 = assemble_seisdata(netstalocchan, fidir, starttime, endtime, data_contents_fraction=0.0)
	println(t1)

	compute_psdpdf(S1, fodir, segments_length, overlap=overlap, smoothing_window_len=smoothing_window_len,
	       figsize=figsize, fmt=fmt, xlims=xlims, ylims=ylims, clims=clims)

end

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
