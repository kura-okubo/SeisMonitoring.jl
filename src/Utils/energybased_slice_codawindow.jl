using SeisNoise, DSP
export energybased_slice_codawindow!, energybased_slice_codawindow

@doc """
    energybased_slice_codawindow!(
			A::AbstractArray, # Time series of correlation function
			maxlag::Real, # maximum time lag [s]
			fm::Real, # mean frequency of time series [Hz]
			fs::Real, # sampling frequency [Hz]
			dist::Real, # distance between station pairs [m]
			background_vel::Real; # background velocity [m/s]
			min_ballistic_twin::Real=1.0, # explicit ballistic time window (see doc)
			nondim_codamaxlag::Real=60, #nondimensional coda lag where kinetic energy is computed.
			nondim_min_coda_length::Real=1.0, #nondimensional minimum coda window length
			nondim_max_coda_length::Real=30.0, #nondimensional maximum coda window length
			coda_energy_threshold::Real=0.95)

Slicing coda window based on coda kinetic energy.

# Arguments
- `A::AbstractArray`: Time series of correlation function
- `maxlag::Real`: maximum time lag [s]
- `fm::Real`: mean frequency of time series [Hz]
- `fs::Real`: sampling frequency [Hz]
- `dist::Real`: distance between station pairs [m]
- `background_vel::Real`: background velocity [m/s]
- `min_ballistic_twin::Real`: explicit ballistic time window (see doc)
- `nondim_codamaxlag::Real=60.0`:nondimensional coda lag where kinetic energy is computed [s]
- 'nondim_min_coda_length::Real=1.0': #nondimensional minimum coda window length
- 'nondim_max_coda_length::Real=30.0': #nondimensional maximum coda window length
- `coda_energy_threshold::Real=0.95`:threshold for the coda endtime based of the fraction of total kinetic energy (see doc)
)
# Return
- `coda_window`: index of coda_window
- `timelag`: time lag [s]
- `fillbox`: coda-window boundaries for plotting coda-window box

# Definision of coda window used in this function
We follow the definision of coda proposed by Aki and Chouet (1975).
Coda window is chosen based on the kinetic energy. We compute I as following:

I(t) = \\int_{0}^{maximum(t)} v(t)^2 dt

This curve looks like shifted sigmoid function. We then normalize I by its maxumum value,
so I(t) ∈ [0, 1.0].

The time for ballistic wave arrival is obtained by max(min_ballistic_twin, tD = r/background_vel), where r is
distance between station. min_ballistic_twin is used to avoid ballistic wave arrival time to be zero in
auto-correlation.

The end of coda window is computed when I(t) reaches coda_energy_threshold (e.g. 0.95), which indicates that the
rest of window is contributed as the background noise.

Hence the time window is sliced by three; 0 < ballistic wave arrival time < end of coda window < maximum(t)

We then evaluate the consistency of these time window as follows:

1. check if the time is ordered as 0 < ballistic wave arrival time < end of coda window < maximum(t).

If not, for example the second and third are flipped, it indicates that there is no coda decay in signal.

2. check the length of coda window such that nondim_min_coda_length< length(coda_window) < nondim_max_coda_length  &&

This ensures that the coda window has enough length, and is not too long to estimate exponential decay.

Author: Kurama Okubo (https://github.com/kura-okubo)
2020.06.25
"""
function energybased_slice_codawindow!(
		A::AbstractArray, # Time series of correlation function
		maxlag::Real, # maximum time lag [s]
		fm::Real, # mean frequency of time series [Hz]
		fs::Real, # sampling frequency [Hz]
		dist::Real, # distance between station pairs [m]
		background_vel::Real; # background velocity [m/s]
		min_ballistic_twin::Real=1.0, # explicit ballistic time window (see doc)
		nondim_codamaxlag::Real=60, #nondimensional minimum coda lag
		nondim_min_coda_length::Real=1.0, #nondimensional minimum coda window length
		nondim_max_coda_length::Real=30.0, #nondimensional maximum coda window length
		coda_energy_threshold::Real=0.95)


	timelag = collect(-maxlag:1/fs:maxlag)

	coda_maxlag = min(maxlag, nondim_codamaxlag/fm)
	mincodalen = nondim_min_coda_length/fm
	maxcodalen = nondim_max_coda_length/fm

	#1. split cc to positive and negative part
	t_pos, t_neg, x_pos, x_neg = split_cc(A, timelag)

	#2. compute minimum coda window at ballistic wave arrival time
	ballistic_arrivaltime = dist / background_vel
	id_tb_pos = findfirst(x-> x >= ballistic_arrivaltime, t_pos)
	id_tb_neg = findfirst(x-> x >= ballistic_arrivaltime, t_neg)

	if any(isnan.([id_tb_pos, id_tb_neg]))
		# discard this cc
		return ([], [], [], Dict())
	end

	t_EB_pos = t_pos[id_tb_pos]
	t_EB_neg = t_neg[id_tb_neg]

	min_coda_pos = max(t_EB_pos, min_ballistic_twin)
	min_coda_neg = max(t_EB_neg, min_ballistic_twin)

	#3. compute cumulative kinetic energy from ballistic wave arrival
	min_coda_pos_id = findfirst(x-> x >= min_coda_pos, t_pos)
	min_coda_neg_id = findfirst(x-> x >= min_coda_neg, t_neg)
	coda_maxlag_pos_id = findfirst(x-> x >= coda_maxlag, t_pos)
	coda_maxlag_neg_id = findfirst(x-> x >= coda_maxlag, t_neg)

	x_pos[1:min_coda_pos_id-1] .= 0.0
	x_neg[1:min_coda_neg_id-1] .= 0.0
	x_pos[coda_maxlag_pos_id:end] .= 0.0
	x_neg[coda_maxlag_neg_id:end] .= 0.0

	#2. compute kinetic energy I(t)
	I_pos = cumtrapz(t_pos, x_pos.^2)
	I_neg = cumtrapz(t_neg, x_neg.^2)

	#3. normalize integrated trace
	I_pos ./= maximum(I_pos)
	I_neg ./= maximum(I_neg)

	#4. compute fraction of ballistic wave arrival and integrated
	id_tc_pos = findfirst(x-> x >= coda_energy_threshold, I_pos)
	id_tc_neg = findfirst(x-> x >= coda_energy_threshold, I_neg)

	if any(isnan.([id_tc_pos, id_tc_neg]))
		# discard this cc
		return ([], [], [], Dict())
	end

	t_EC_pos = t_pos[id_tc_pos]
	t_EC_neg = t_neg[id_tc_neg]

	EB_pos = I_pos[id_tb_pos]
	EB_neg = I_neg[id_tb_neg]
	EC_pos = I_pos[id_tc_pos]
	EC_neg = I_neg[id_tc_neg]

	max_coda_pos = t_EC_pos
	max_coda_neg = t_EC_neg
	# sanity check of coda window
	coda_window_len_pos = max_coda_pos - min_coda_pos
	coda_window_len_neg = max_coda_neg - min_coda_neg

	#1. order of time window

	flag_coda_pos = (0 < min_coda_pos < max_coda_pos) && # check if order of coda window is consistent
	 				(mincodalen <= coda_window_len_pos <= maxcodalen) ? true : false # check if end of coda is less than coda_limtaion_α*min_coda_pos

	flag_coda_neg = (0 < min_coda_neg < max_coda_neg) &&
	 				(mincodalen <= coda_window_len_neg <= maxcodalen) ? true : false


	# compute coda window id and fillbox for plotting
	coda_window = Int[]
	flag_coda_pos && (coda_window = vcat(coda_window, findall(x -> min_coda_pos<= x <= max_coda_pos, timelag)))
	flag_coda_neg && (coda_window = vcat(coda_window, findall(x -> -max_coda_neg<= x <= -min_coda_neg, timelag))) # note for the sign of time

	fill_box = []
	flag_coda_neg && (fill_box = vcat(fill_box, [-max_coda_neg, -min_coda_neg]))
	flag_coda_pos && (fill_box = vcat(fill_box, [min_coda_pos, max_coda_pos]))

	CodaSliceDict=Dict(
		"t_pos" => t_pos,
		"t_neg" => t_neg,
		"x_pos" => x_pos,
		"x_neg" => x_neg,
		"I_pos" => I_pos,
		"I_neg" => I_neg,
		"t_EB_pos" => min_coda_pos,
		"t_EB_neg" => min_coda_neg,
		"t_EC_pos" => max_coda_pos,
		"t_EC_neg" => max_coda_neg,
		"EB_pos" => EB_pos,
		"EB_neg" => EB_neg,
		"EC_pos" => EC_pos,
		"EC_neg" => EC_neg,)

	#evaluate coda window and return values
	if isempty(coda_window)
		@warn("Coda window is null. Please check parameters.")
		return ([], [], [], CodaSliceDict)
	end

	return (coda_window, timelag, fill_box, CodaSliceDict)

end


"""
	cumtrapz(X::AbstractVector, Y::AbstractVector)

Compute accumulate integration of X.
"""
function cumtrapz(X::AbstractVector, Y::AbstractVector)
  # Check matching vector length
  @assert length(X) == length(Y)
  # Initialize Output
  out = similar(X)
  out[1] = 0
  # Iterate over arrays
  for i in 2:length(X)
    out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
  end
  # Return output
  return out
end


@doc (@doc slice_codawindow!)
function energybased_slice_codawindow(
	A::AbstractArray, # Time series of correlation function
	maxlag::Real, # maximum time lag [s]
	fm::Real, # mean frequency of time series [Hz]
	fs::Real, # sampling frequency [Hz]
	dist::Real, # distance between station pairs [m]
	background_vel::Real; # background velocity [m/s]
	min_ballistic_twin::Real=1.0, # explicit ballistic time window (see doc)
	nondim_codamaxlag::Real=60, #nondimensional minimum coda lag
	nondim_min_coda_length::Real=1.0, #nondimensional minimum coda window length
	nondim_max_coda_length::Real=30.0, #nondimensional maximum coda window length
	coda_energy_threshold::Real=0.95)

	U = deepcopy(A);
	coda_window, timelag, fillbox, CodaSliceDict = energybased_slice_codawindow!(U,maxlag,fm,fs,dist,
	background_vel,
	min_ballistic_twin=min_ballistic_twin,
	nondim_codamaxlag=nondim_codamaxlag,
	nondim_min_coda_length=nondim_min_coda_length,
	nondim_max_coda_length=nondim_max_coda_length,
	coda_energy_threshold=coda_energy_threshold)

	return (U, coda_window, timelag, fillbox, CodaSliceDict)
end

@doc (@doc slice_codawindow!)
function energybased_slice_codawindow!(
	C::CorrData,
	background_vel::Real; # background velocity [m/s]
	min_ballistic_twin::Real=3.0, # explicit ballistic time window (see doc)
	nondim_codamaxlag::Real=60, #nondimensional minimum coda lag
	nondim_min_coda_length::Real=5.0,
	nondim_max_coda_length::Real=3.0,
	coda_energy_threshold::Real=0.95)

	# C.dist is in km, thus multiplied by 1e3
	coda_window, timelag, fillbox, CodaSliceDict = energybased_slice_codawindow!(
	C.corr,C.maxlag,(C.freqmax+C.freqmin)/2.0, C.fs,C.dist*1e3,background_vel,
	min_ballistic_twin=min_ballistic_twin,
	nondim_codamaxlag=nondim_codamaxlag, #nondimensional minimum coda lag
	nondim_min_coda_length=nondim_min_coda_length,
	nondim_max_coda_length=nondim_max_coda_length,
	coda_energy_threshold=coda_energy_threshold)

	return (coda_window, timelag, fillbox, CodaSliceDict)
end
