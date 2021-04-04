using SeisMonitoring: assemble_corrdata, cc_medianmute!, smooth_withfiltfilt
using SeisIO, SeisNoise, JLD2, Dates, Plots, Statistics, DSP, Interpolations

@doc """
    const_slice_codawindow!(A, maxlag, fm, fs, dist, background_vel,
	min_ballistic_twin, max_coda_length; zeropad=false)

Slicing coda window based on prescribed velocity, maximum coda length and maximum coda length.

Author: Kurama Okubo (https://github.com/kura-okubo)
2021.03.15
"""
function const_slice_codawindow(
		A::AbstractArray, 			# set of time series of correlation function
		maxlag::Real, 				# maximum time lag [s]
		fm::Real, 					# mean frequency of time series [Hz]
		fs::Real, 					# sampling frequency [Hz]
		dist::Real, 				# distance between station pairs [m]
		background_vel::Real, 		# background velocity [m/s]
		min_ballistic_twin::Real, 	# explicit ballistic time window (see doc)
		max_coda_length::Real; 		# maximum coda window length [s]
		coda_init_factor::Real=2.0, # coda window starts from coda_init_factor*dist/vel
		coda_minlen_factor::Real=5.0 # minimumlength is determined by this factor * (1/fm, period of cc) * fs points.
		)

	# central index (zero time lag)
	timelag = collect(-maxlag:1/fs:maxlag)
	Ntimelag = length(timelag)

	if length(timelag) != size(A, 1)
		@warn("Length of time series is not consitent with max time lag. skipping.")
		return ([], [], [], Dict())
	end

	centerid = round(Int, Ntimelag/2)

	#1. Minimum explicit ballistic time window
	min_ballistic_width = min_ballistic_twin * fs
	minbal_window_neg =  round(Int, centerid - min_ballistic_width)
	minbal_window_pos =  round(Int, centerid + min_ballistic_width)
	minbal_window = vcat(collect(1:minbal_window_neg), collect(minbal_window_pos:Ntimelag))

	#2. Distance dependent window and max window
	if round(Int, fs * dist / background_vel) < 1
		# this is the case for auto-correlation
		# Ballistic time window approximated by background velocity is entire time lag in this case.
		ba_window_neg = centerid
		ba_window_pos = centerid
		# Maximum coda window
		max_coda_window_neg = round(Int, centerid - max_coda_length * fs)
		max_coda_window_pos = round(Int, centerid + max_coda_length * fs)

		# align coda window border if it is outside of tvec
		max_coda_window_neg < 1 	   && (max_coda_window_neg = 1)
		max_coda_window_pos > Ntimelag && (max_coda_window_pos = Ntimelag)

		approx_coda_window = collect(max_coda_window_neg:max_coda_window_pos)

	else
	   # this is the case for cross-correlation between station pair
	   # Ballistic time window approximated by background velocity is removed from the coda window in this case.
	   ba_window_neg = round(Int, centerid - coda_init_factor * fs * dist / background_vel) #[m] / [m/s]
	   ba_window_pos = round(Int, centerid + coda_init_factor * fs * dist / background_vel) #[m] / [m/s]

	   # Maximum coda window
	   max_coda_window_neg =  round(Int, centerid - max_coda_length * fs) #[m] / [m/s]
	   max_coda_window_pos =  round(Int, centerid + max_coda_length * fs) #[m] / [m/s]

	   # align coda window border if it is outside of tvec
	   ba_window_neg < 1 		&& (ba_window_neg = 1)
	   ba_window_pos > Ntimelag && (ba_window_pos = Ntimelag)
	   max_coda_window_neg < 1 	&& (max_coda_window_neg = 1)
	   max_coda_window_pos > Ntimelag && (max_coda_window_pos = Ntimelag)

	   coda_window_neg = collect(max_coda_window_neg:ba_window_neg)
	   coda_window_pos = collect(ba_window_pos:max_coda_window_pos)

	   approx_coda_window = vcat(coda_window_neg, coda_window_pos)

	end

	coda_window = intersect(minbal_window, approx_coda_window)

	CodaSliceDict=Dict("timelag" => timelag)

	# Evaluate if positive and negative coda_window length is longer than minimum threshold
	# 1. decompose coda window into positive and negative
	coda_pos_ind = findall(x -> x>=centerid , coda_window)
	coda_neg_ind = findall(x -> x<centerid , coda_window)
	# 2. compoute minimum length of one-side coda window
	min_codalength_points = trunc(Int, coda_minlen_factor * fs / fm)

	# Check if coda window is symmetric
	len_coda_pos = length(coda_pos_ind)
	len_coda_neg = length(coda_neg_ind)
	if (len_coda_pos != length(coda_neg_ind))
		@warn("Coda window is not symmetric. Num. of coda Pos:$(coda_pos_ind), Neg:$(len_coda_neg)")
	else
		len_coda_window = len_coda_pos
	end

	# Check if coda window has enought length
	if ( len_coda_window < min_codalength_points )
		# both positive and negative do not have enough coda window length
		@warn("Coda window does not have enought length. Skippling")
		return ([], [], [], CodaSliceDict)
	else
		# both positive and negative have enough coda window length
		coda_window_all = coda_window
	end

	# NOTE: fill_box will be deprecated due to the change of coda definision.
	fill_box = [timelag[minimum(coda_window_all)], timelag[min(ba_window_neg,minbal_window_neg)],
				timelag[max(ba_window_pos,minbal_window_pos)], timelag[maximum(coda_window_all)]]

	return (coda_window_all, timelag, fill_box, CodaSliceDict)
end

@doc (@doc const_slice_codawindow)
function const_slice_codawindow(
		C::CorrData, # set of time series of correlation function
		background_vel::Real, # background velocity [m/s]
		min_ballistic_twin::Real, # explicit ballistic time window (see doc)
		max_coda_length::Real; # maximum coda window length [s]
		coda_init_factor::Real=3.0, # coda window starts from coda_init_factor*dist/vel
		coda_minlen_factor::Real=5.0, # minimum coda length defined as min_codalength_α * mwcc_len
		debugplot::Bool=false, # plot debug figures
		foname::String="", # figure name for debug plot
		fodir::String="",
		xlims::Tuple=(-40, 40)
		)

	# C.dist is in km, thus multiplied by 1e3
	coda_window, timelag, fillbox, CodaSliceDict = const_slice_codawindow(
	C.corr, C.maxlag, (C.freqmax+C.freqmin)/2.0, C.fs, C.dist*1e3, background_vel,
	min_ballistic_twin, max_coda_length, coda_init_factor=coda_init_factor,
	coda_minlen_factor=coda_minlen_factor)

	# debug plot is available only with corr data.
	if debugplot
		p1 = corrplot(C)
		xlims!(xlims)
		# plot vlines on coda window
		!isempty(coda_window) && (p1 = vline!(timelag[coda_window], width=1.0, color=:orange, legend=false, alpha=0.3))

		xlims!(xlims)
		xlabel!("Time lag[s]")
		ylabel!("coherence")
		title!("$(C.name) $(C.id)")

		plot(p1, size=(1200,400))

		fopath = joinpath(fodir, foname)*".png"
		savefig(fopath)
	end

	return (coda_window, timelag, fillbox, CodaSliceDict)
end
#
# #
# # Test script
# #1. Cross-correlation
# finame = "./cc_data/BP.SCYB.40.SP1-BP.SMNB.40.SP2__2013-10-25T00:00:00__2013-10-30T00:00:00.jld2"
# fi = jldopen(finame)
# CorrData_Buffer = Dict()
# starttime = DateTime(2013,10,25)
# endtime = DateTime(2013,10,26)
# freqkey = "0.9-1.2"
# C, CorrData_Buffer = assemble_corrdata(fi,starttime,endtime,freqkey,
#                                         min_cc_datafraction = 0.5)
# cc_medianmute!(C, 2.0, 0.1)
#
# background_vel=1000
# min_ballistic_twin=1.0
# coda_init_factor=2.0
# max_coda_length=40.0
# coda_minlen_factor = 5.0
#
# stationpairname = splitdir(finame)[2][1:end-5]
# foname = stationpairname*"_"*freqkey
# fodir = "./debugplot"
# !ispath(fodir) && mkpath(fodir)
# xlims = (-80, 80)
#
# coda_window, timelag, fillbox, CodaSliceDict = const_slice_codawindow(C,background_vel, min_ballistic_twin, max_coda_length,
# 		coda_init_factor=coda_init_factor, coda_minlen_factor=coda_minlen_factor,
# 		debugplot=true, foname=foname, fodir=fodir, xlims=xlims)
# close(fi)
# #2. Auto-correlation
# finame = "./cc_data/BP.SCYB.40.SP1-BP.SCYB.40.SP1__2013-10-25T00:00:00__2013-10-30T00:00:00.jld2"
# fi = jldopen(finame)
# CorrData_Buffer = Dict()
# starttime = DateTime(2013,10,25)
# endtime = DateTime(2013,10,26)
# freqkey = "0.9-1.2"
# C, CorrData_Buffer = assemble_corrdata(fi,starttime,endtime,freqkey,
#                                         min_cc_datafraction = 0.5)
# cc_medianmute!(C, 2.0, 0.1)
#
# background_vel=1000
# min_ballistic_twin=3.0
# coda_init_factor=2.0
# max_coda_length=40.0
# coda_minlen_factor = 5.0
#
# stationpairname = splitdir(finame)[2][1:end-5]
# foname = stationpairname*"_"*freqkey
# fodir = "./debugplot"
# !ispath(fodir) && mkpath(fodir)
# xlims = (-80, 80)
#
# coda_window, timelag, fillbox, CodaSliceDict = const_slice_codawindow(C,background_vel, min_ballistic_twin, max_coda_length,
# 		coda_init_factor=coda_init_factor, coda_minlen_factor=coda_minlen_factor,
# 		debugplot=true, foname=foname, fodir=fodir, xlims=xlims)
