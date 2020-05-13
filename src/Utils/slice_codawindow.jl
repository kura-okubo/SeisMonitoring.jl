using SeisNoise, DSP
export slice_codawindow!, slice_codawindow
@doc """
    slice_codawindow!(A, maxlag, fm, fs, dist, background_vel, coda_Qinv,
	min_ballistic_twin, max_coda_length; attenuation_minthreshold=0.1, zeropad=false)

Slicing coda window based on prescribed velocity and coda-Q^{-1}.

# Arguments
- `A::AbstractArray`: Time series of correlation function
- `maxlag::Real`: maximum time lag [s]
- `fm::Real`: mean frequency of time series [Hz]
- `fs::Real`: sampling frequency [Hz]
- `dist::Real`: distance between station pairs [m]
- `background_vel::Real`: background velocity [m/s]
- `coda_Qinv::Real`: coda Q inverse; set 0 if not using attenuation threshold
- `min_ballistic_twin::Real`: explicit ballistic time window (see doc)
- `max_coda_length::Real`:maximum coda window length [s]
- `attenuation_minthreshold::Real=0.1`:threshold for attenuation decay (see doc)
- `zeropad::Bool=false`: zero padding outside of coda window
)
# Return
- `coda_window`: index of coda_window
- `timelag`: time lag [s]
- `fillbox`: coda-window boundaries for plotting coda-window box

# Definision of coda window used in this function
We follow the definision of coda proposed by Aki and Chouet (1975).
The local energy in coda part comes from various azimuths, showing that
the origin of coda wave is not source, but scattered wave due to
structural heterogeneity of material.

It is hard, and not essential, to find out the border between ballistic (i.e.
direct arrival) and coda wave due to uncertanity of ballistic wave arrival.
In addition, it might be overlapped with given source spectrum.
Thus, we generally approximate the border with a factor
of station distance devided by background velocity.

It is noteworthy that the detailed background velocity and coda_Qinv
are not required in essence; the representative values are sufficient
to approximate coda window.

For the autocorrelation, entire timelag is coda in principal. However,
we avoid around zero timelag by explicit minimum ballistic window length.

#The schematic of coda window slicing performed in this function

Given positive correlation function as below:

```
|     .-.
|    /   \\         .-.
|   /     \\       /   \\       .-.     .-.     _   _
+--/-------\\-----/-----\\-----/---\\---/---\\---/-\\-/-\\/\\/--->
| /         \\   /       \\   /     '-'     '-'
|/           '-'         '-'
0
Time lag [s]
```

1. Minimum explicit ballistic time window
```
|         |-----------------------------------------
0	min_ballistic_timewindow[s]
```

2. coda_Q attenuation
*Note that we use normalized coda_Q exponential decay.
Not performing logscale fitting on the data.

```
1.0 +- --                               |
    |     ---                           |
    |        -----                      |
    |             ------                |
α   +-                  ----------      |
    |                             ------|
|---------------------------------|
0                          attenuation_timewindow[s]
```

- `α`: attenuation_minthreshold (default: 0.1)
- `normalized coda_Q_decaymodel(t) = exp.((-pi*fm*coda_Qinv) .* abs(t))``
- `attenuation_timewindow = filter(t -> coda_Q_decaymodel(t) > α, timelag)`

3. Ballistic time window approximated by background velocity
```
|               |-----------------------------------
0	ballistic_arrivaltime_factor * (station distance/background_vel[m/s])
```

4. Maximum coda window
```
|------------------------------------------|
0						      max_coda_timewindow[s]
```

The final coda window is computed as an intersect of the four criteria above:
```
|               |=================|
0        		  pos_coda_window
```

Author: Kurama Okubo (https://github.com/kura-okubo)
2020.04.29
"""
function slice_codawindow!(
		A::AbstractArray, # Time series of correlation function
		maxlag::Real, # maximum time lag [s]
		fm::Real, # mean frequency of time series [Hz]
		fs::Real, # sampling frequency [Hz]
		dist::Real, # distance between station pairs [m]
		background_vel::Real, # background velocity [m/s]
		coda_Qinv::Real, # coda Q inverse; set 0 if not using attenuation threshold
		min_ballistic_twin::Real, # explicit ballistic time window (see doc)
		max_coda_length::Real; # maximum coda window length [s]
		attenuation_minthreshold::Real=0.1,
		zeropad::Bool=false
		)

	# central index (zero time lag)
	timelag = collect(-maxlag:1/fs:maxlag)
	Ntimelag = length(timelag)

	if length(timelag) != size(A, 1)
		@warn("Length of time series is not consitent with max time lag. skipping.")
		return (NaN, NaN, NaN)
	end

	centerid = round(Int, Ntimelag/2)

	#1. Minimum explicit ballistic time window
	min_ballistic_width = min_ballistic_twin * fs
	minbal_window_neg =  round(Int, centerid - min_ballistic_width)
	minbal_window_pos =  round(Int, centerid + min_ballistic_width)
	minbal_window = vcat(collect(1:minbal_window_neg), collect(minbal_window_pos:Ntimelag))

	#2. coda_Q attenuation
	coda_Q_decaymodel(t) = exp((-pi*fm*coda_Qinv) * abs(t))
	attenuation_timewindow = findall(t -> coda_Q_decaymodel(t) > attenuation_minthreshold, timelag)

	if round(Int, fs * dist / background_vel) < 1
		# this is the case for auto-correlation
		#3. Ballistic time window approximated by background velocity is entire time lag in this case.
		ba_window_neg = centerid
		ba_window_pos = centerid
		#4. Maximum coda window
		max_coda_window_neg = round(Int, centerid - max_coda_length * fs)
		max_coda_window_pos = round(Int, centerid + max_coda_length * fs)

		# align coda window border if it is outside of tvec
		max_coda_window_neg < 1 	   && (max_coda_window_neg = 1)
		max_coda_window_pos > Ntimelag && (max_coda_window_pos = Ntimelag)

		approx_coda_window = collect(max_coda_window_neg:max_coda_window_pos)

	else
	   # this is the case for cross-correlation between station pair
	   #3. Ballistic time window approximated by background velocity is null in this case.
	   ba_window_neg = round(Int, centerid - fs * dist / background_vel) #[m] / [m/s]
	   ba_window_pos = round(Int, centerid + fs * dist / background_vel) #[m] / [m/s]

	   #4. Maximum coda window
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

	coda_window = intersect(minbal_window,
							attenuation_timewindow,
							approx_coda_window)

	#evaluate coda window and return values
	if isempty(coda_window)
	    @warn("Coda window is null. Please check parameters.")
		return (NaN, NaN, NaN)
	end

	# zero padding outside of coda window if zeropad == true
	if zeropad
		coda_neg_id = minimum(coda_window):min(ba_window_neg,minbal_window_neg)
		coda_pos_id = max(ba_window_pos,minbal_window_pos):maximum(coda_window)
		# apply tukey window within the coda window
		tneg = DSP.tukey(length(coda_neg_id), 0.1)
		tpos = DSP.tukey(length(coda_pos_id), 0.1)
		for i = 1:size(A, 2)
			A[filter(x -> !(x in coda_window), 1:Ntimelag), i] .= 0.0
			A[coda_neg_id,i] .*= tneg
			A[coda_pos_id,i] .*= tpos
		end
	end

	fill_box = [timelag[minimum(coda_window)], timelag[min(ba_window_neg,minbal_window_neg)],
				timelag[max(ba_window_pos,minbal_window_pos)], timelag[maximum(coda_window)]]

	return (coda_window, timelag, fill_box)
end

@doc (@doc slice_codawindow!)
function slice_codawindow(
	A::AbstractArray, # Time series of correlation function
	maxlag::Real, # maximum time lag [s]
	fm::Real, # mean frequency of time series [fs]
	fs::Real, # sampling frequency [fs]
	dist::Real, # distance between station pairs [m]
	background_vel::Real, # background velocity [m/s]
	coda_Qinv::Real, # coda Q inverse; set 0 if not using attenuation threshold
	min_ballistic_twin::Real, # explicit ballistic time window (see doc)
	max_coda_length::Real; # maximum coda window length [s]
	attenuation_minthreshold::Real=0.1,
	zeropad::Bool=false
	)
	U = deepcopy(A);
	coda_window, timelag, fillbox = slice_codawindow!(U,maxlag,fm,fs,dist,
	background_vel,coda_Qinv,min_ballistic_twin,max_coda_length,
	attenuation_minthreshold=attenuation_minthreshold,
	zeropad=zeropad)
	return (U, coda_window, timelag, fillbox)
end

@doc (@doc slice_codawindow!)
function slice_codawindow!(
	C::CorrData,
	background_vel::Real, # background velocity [m/s]
	coda_Qinv::Real, # coda Q inverse; set 0 if not using attenuation threshold
	min_ballistic_twin::Real, # explicit ballistic time window (see doc)
	max_coda_length::Real; # maximum coda window length [s]
	attenuation_minthreshold::Real=0.1,
	zeropad::Bool=false
	)
	# C.dist is in km, thus multiplied by 1e3
	coda_window, timelag, fillbox = slice_codawindow!(
	C.corr,C.maxlag,(C.freqmax+C.freqmin)/2.0, C.fs,C.dist*1e3,background_vel,coda_Qinv,
	min_ballistic_twin,max_coda_length,
	attenuation_minthreshold=attenuation_minthreshold, zeropad=zeropad);
	return (coda_window, timelag, fillbox)
end
