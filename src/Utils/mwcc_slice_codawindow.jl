using SeisMonitoring: assemble_corrdata, cc_medianmute!, smooth_withfiltfilt
using SeisIO, SeisNoise, JLD2, Dates, Plots, Statistics, DSP, Interpolations

export mwcc_slice_codawindow

@doc """
    mwcc_slice_codawindow!(A, maxlag, fm, fs, dist, background_vel, coda_Qinv,
	min_ballistic_twin, max_coda_length; attenuation_minthreshold=0.1, zeropad=false)

Slicing coda window based on prescribed velocity and correlation.

Author: Kurama Okubo (https://github.com/kura-okubo)
2020.07.24
"""
function mwcc_slice_codawindow(
		A::AbstractArray, # set of time series of correlation function
		maxlag::Real, # maximum time lag [s]
		fm::Real, # mean frequency of time series [Hz]
		fs::Real, # sampling frequency [Hz]
		dist::Real, # distance between station pairs [m]
		background_vel::Real, # background velocity [m/s]
		min_ballistic_twin::Real, # explicit ballistic time window (see doc)
		max_coda_length::Real; # maximum coda window length [s]
		mwcc_threshold::Real=0.5, # threshold of mwcc
		coda_init_factor::Real=2.0, # coda window starts from coda_init_factor*dist/vel
		mwcc_len_α::Real=3.0, # mwcc length is defined as mwcc_len_α*fm
		# zeropad::Bool=false
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
	   # Ballistic time window approximated by background velocity is null in this case.
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

	#3. compute moving window correlation coefficient
	mwcc_len = trunc(Int, (mwcc_len_α/fm)*fs)
	mwcc_step = trunc(Int, mwcc_len/2)
	cc = compute_meanmwcc(A, mwcc_len, mwcc_step)
	# apply smoothing
	cc = smooth_withfiltfilt(cc, window_len=5*mwcc_len)

	mwcc_window = findall(x -> x > mwcc_threshold, cc)

	coda_window = intersect(minbal_window,
							approx_coda_window,
							mwcc_window)

	CodaSliceDict=Dict("cc" => cc,
						"timelag" => timelag)

	#evaluate coda window and return values
	if isempty(coda_window)
	    @warn("Coda window is null. Please check parameters.")
		return ([], [], [], CodaSliceDict)
	end

	# NOTE: fill_box will be deprecated due to the change of coda definision.
	fill_box = [timelag[minimum(coda_window)], timelag[min(ba_window_neg,minbal_window_neg)],
				timelag[max(ba_window_pos,minbal_window_pos)], timelag[maximum(coda_window)]]

	return (coda_window, timelag, fill_box, CodaSliceDict)
end


"""
    mwcc(x::Vector, y::Vector, len::Int, step::Int)

moving window correlation coefficient
# Input
- x, y: two traces with same length
- len: length of moving window (points)
- step: step of moving window (points)

# Output
- cc: vector of correlation coefficients interpolated to be the same length as x and y

# Tips

1. len can be set based on factor of mean frequency of x and y. e.g. len = trunc(Int, 3 * fm).

2. step can be half of len. e.g. step = trunc(Int, len/2)

3. To keep symmetricity, compute cc from both sides of data
"""
function mwcc(x::AbstractArray, y::AbstractArray, len::Int, step::Int)
    length(x) != length(y) && error("mwcc: x and y should be same length")
    N = length(x)
    st = range(1, stop=N-len, step=step)
    mt = convert(Int, round(len/2))
    tt = range(mt, stop=maximum(st)+mt, step=step)

	# forward cc
	ct = []
    for is in st
        es=is+len
        xt = @view x[is:es]
        yt = @view y[is:es]
        push!(ct, cor(xt, yt))
    end
    # linear interpolate to get N length correlation coefficient
    itp = Interpolations.interpolate(ct, BSpline(Linear()))
    etpf = Interpolations.extrapolate(itp, Flat())
    sitp = Interpolations.scale(etpf, tt)
    cc = sitp(1:N)

	# inverted cc
	x_inv = x[end:-1:1]
	y_inv = y[end:-1:1]
	ct_inv = []
	for is in st
		es=is+len
		xt = @view x_inv[is:es]
		yt = @view y_inv[is:es]
		push!(ct_inv, cor(xt, yt))
	end
	# linear interpolate to get N length correlation coefficient
	itp_inv = Interpolations.interpolate(ct_inv, BSpline(Linear()))
	etpf_inv = Interpolations.extrapolate(itp_inv, Flat())
	sitp_inv = Interpolations.scale(etpf_inv, tt)
	cc_inv = sitp(1:N)

    return (cc .+ cc_inv[end:-1:1])./2
end

function compute_meanmwcc(A::AbstractArray, len::Int, step::Int)

	# compute pairwise moving window correlation
	T, N = size(A)
	Npairwise = trunc(Int, N*(N-1)/2)
	CC_sum = zeros(Float32, T)

	@inbounds for i = 2:N
	    for j = 1:i-1
	        x = @view A[:, i]
	        y = @view A[:, j]
	        cc = mwcc(x, y, len, step)
	        # avoid nan value
	        replace!(cc, NaN=>-1)
	        # any(isnan.(cc)) && continue
	        CC_sum += cc
	    end
	end

	return CC_sum ./ Npairwise
end

@doc (@doc mwcc_slice_codawindow)
function mwcc_slice_codawindow(
		C::CorrData, # set of time series of correlation function
		background_vel::Real, # background velocity [m/s]
		min_ballistic_twin::Real, # explicit ballistic time window (see doc)
		max_coda_length::Real; # maximum coda window length [s]
		mwcc_threshold::Real=0.5, # threshold of mwcc
		coda_init_factor::Real=2.0, # coda window starts from coda_init_factor*dist/vel
		mwcc_len_α::Real=3.0, # mwcc length is defined as mwcc_len_α*fm
		debugplot::Bool=false, # plot debug figures
		foname::String="", # figure name for debug plot
		fodir::String="",
		xlims::Tuple=(-40, 40)
		# zeropad::Bool=false
		)

	# C.dist is in km, thus multiplied by 1e3
	coda_window, timelag, fillbox, CodaSliceDict = mwcc_slice_codawindow(
	C.corr, C.maxlag, (C.freqmax+C.freqmin)/2.0, C.fs, C.dist*1e3, background_vel,
	min_ballistic_twin, max_coda_length, mwcc_threshold=mwcc_threshold,
	coda_init_factor=coda_init_factor, mwcc_len_α=mwcc_len_α)

	# debug plot is available only with corr data.
	if debugplot && !isempty(coda_window)
		p1 = corrplot(C)
		xlims!(xlims)
		# plot vlines on coda window
		p1 = vline!(timelag[coda_window], width=1.0, color=:orange, legend=false, alpha=0.3)

		p2 = plot(timelag, CodaSliceDict["cc"], label="")
		p2 = vline!(timelag[coda_window], width=1.0, color=:orange, legend=false, alpha=0.3)

		xlims!(xlims)
		xlabel!("Time lag[s]")
		ylabel!("cc")
		title!("Moving window correlation coefficient")

		plot(p1, p2, layout=grid(2, 1, heights=[0.7, 0.3]), size=(1200,800))

		fopath = joinpath(fodir, foname)*".png"
		savefig(fopath)
	end

	return (coda_window, timelag, fillbox, CodaSliceDict)
end
#
# # Test script
# station = "CLC" #CGO
# fi = jldopen("./CI.$(station)-CI.$(station).jld2")
# CorrData_Buffer = Dict()
# stachanpair = "CI.$(station)..BHZ-CI.$(station)..BHZ"
# starttime = DateTime(2019,4,1)
# endtime = DateTime(2020,4,1)
# freqkey = "1.0-2.0"
# C, CorrData_Buffer = assemble_corrdata(fi,stachanpair,starttime,endtime,freqkey,
#                                         min_cc_datafraction = 0.5)
#
# cc_medianmute!(C, 2.0)
#
# C.dist = 0.0
# background_vel=1000
# min_ballistic_twin=1.0
# max_coda_length=60.0
# mwcc_threshold=0.4
# mwcc_len_α = 3.0
#
# foname = C.name
# fodir = "./debugplot"
# !ispath(fodir) && mkpath(fodir)
#
# coda_window, timelag, fillbox, CodaSliceDict = mwcc_slice_codawindow(C,background_vel, min_ballistic_twin, max_coda_length,
# 		mwcc_threshold=mwcc_threshold, coda_init_factor=2.0, mwcc_len_α=mwcc_len_α, debugplot=true, foname=foname, fodir=fodir)
