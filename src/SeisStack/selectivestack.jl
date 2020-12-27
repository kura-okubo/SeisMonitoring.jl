"""
    selectivestack!(C::CorrData, stackmode::String; dist_threshold::AbstractFloat=0.0,
    					distance_type::DataType=CorrDist())

Stack the windows in a CorrData object that exceed a correlation-coefficient threshold with respect to a reference.

# Arguments
- `C::CorrData` : Input CorrData matrix used to define the reference and each window
- `stackmode::String` : If reference, applying linear stack, while applying selective stack when stack mode is shorttime.
- `dist_threshold::AbstractFloat` : distance threshold for selective stack.
- `distance_type::DataType` : Method to measure distance with respect to a reference. The choice can be found at https://github.com/JuliaStats/Distances.jl.

# Note
If no traces are accepted during seleoction return CorrData with C.corr[:,1] = [NaN, ..., NaN].
"""
function selectivestack!(C::CorrData, stackmode::String; dist_threshold::AbstractFloat=0.0,
					distance_type::String="CorrDist")

	disttype = eval(Meta.parse(distance_type)) 	# parse distance_type to compute distance between reference and current trace.

	if stackmode == "reference"
		# in reference mode, apply linear stack as longterm stack is assumed to be converged
		stack!(C, allstack=true, stacktype=mean)
		return

	elseif stackmode == "shorttime"
		if !haskey(C.misc, "reference")
			@warn("$(C.name) does not have reference curve at selective stack. Thus apply linearstack.")
			stack!(C, allstack=true, stacktype=mean)
			return
		end
		# compute Pearson's correlation coefficient
		T, N = size(C.corr)
		refmatrix = Array{Float32, 2}(undef, T, N) # make refmatrix to use colwise()
		# println(C.misc["reference"])
		for i=1:N
			# DEBUG:
			try
				refmatrix[:, i] = C.misc["reference"]
			catch y
				println(C)
				println(C.misc["reference"])
				println(y)
			end
		end
		r = Distances.colwise(disttype(), C.corr, refmatrix)
		# inds = findall(x -> x > dist_threshold, r) # find all traces which are closer than dist_threshold
		#NOTE: CorrDist is bounded from [0, 2] in Distances.jl
		inds = findall(x -> x < dist_threshold, r)
		C.misc["selectivestack_acceptratio"] = length(inds)/N
		C.misc["distance_type"] = distance_type
		C.corr = C.corr[:,inds]
		# println(size(C.corr))
		stack!(C, allstack=true, stacktype=mean)
		return
	end
end

selectivestack(C::CorrData, stackmode::String; dist_threshold::AbstractFloat=0.0,
distance_type=CorrDist()) = (U = deepcopy(C);
      	 selectivestack!(U,stackmode,dist_threshold=dist_threshold,
             distance_type=distance_type);return U)

#
# # Test
# using SeisIO, SeisNoise, DSP, Statistics, Distances
#
# c1 = jldopen("corr_test.jld2", "r")
# C = c1["test"]
# close(c1)
# # manipulate time
# starttime = DateTime("2004-04-01T00:00:00")
# endtime = DateTime("2004-04-02T00:00:00")
#
# C.misc["reference"] = stack(C, allstack=true).corr
#
# C2 = selectivestack(C, "reference", dist_threshold=0.0,distance_type="Euclidean")
#
# C3 = selectivestack(C, "shorttime", dist_threshold=0.3,distance_type="CorrDist")
#
# selectivestack!(C, "shorttime", dist_threshold=0.1,distance_type="CorrDist")
