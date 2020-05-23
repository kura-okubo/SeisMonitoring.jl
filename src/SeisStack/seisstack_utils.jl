function get_shorttime_window(starttime, endtime, cc_time_unit, averagestack_factor::Int=1, averagestack_step::Int=1)
    su = d2u(starttime)
    eu = d2u(endtime)

	shorttime_window_length = cc_time_unit*averagestack_factor
	shorttime_window_step = cc_time_unit*averagestack_step
    starts = Array(range(su,stop=eu-shorttime_window_length,step=shorttime_window_step))
	ends = starts .+ shorttime_window_length
    return (u2d.(starts), u2d.(ends))
end

# starttime = DateTime("2004-04-01T00:00:00")
# endtime = DateTime("2004-04-03T00:00:00")
# cc_time_unit = 86400
# averagestack_factor = 1
#
# starts, ends = get_shorttime_window(starttime, endtime, cc_time_unit, averagestack_factor)

"""
	get_reference(fi)
return CorrData dictionary of reference traces.
"""
function get_reference(fi)
    ReferenceDict = Dict()
    stachanpairs =  keys(fi)
    for stachanpair in stachanpairs
        reftimekey = keys(fi[stachanpair])
        if !isone(length(reftimekey))
            # reference jld2 file has multiple reference.
            error("reference $(stachanpair) has multiple reference traces. Please check reference jld2 file.")
        end
        for freqkey in keys(fi[joinpath(stachanpair, reftimekey[1])])
            # remove information of reference time window
            refdictpath = joinpath(stachanpair, freqkey) #BP.CCRB..BP1-BP.EADB..BP1/0.1-0.2
            ReferenceDict[refdictpath] = fi[joinpath(stachanpair, reftimekey[1], freqkey)]
        end
    end
    return ReferenceDict
end


"""
    append_reference!(C::CorrData, stachanpair::String, ReferenceDict::Dict)

Append reference trace in C.misc["reference"] from ReferenceDict
"""
function append_reference!(C::CorrData, stachanpair::String, freqkey::String, ReferenceDict::Dict, InputDict::OrderedDict)
    refdictpath = joinpath(stachanpair, freqkey) #BP.CCRB..BP1-BP.EADB..BP1/0.1-0.2
	# println(keys(ReferenceDict))
	# println("debug: refdictpath: $(refdictpath)")

    if haskey(ReferenceDict, refdictpath)
		Ctemp = ReferenceDict[refdictpath]

		slice_codawindow!(Ctemp,
							InputDict["background_vel"],
							InputDict["coda_Qinv"],
							InputDict["min_ballistic_twin"],
							InputDict["max_coda_length"],
							attenuation_minthreshold=InputDict["slice_minthreshold"],
							zeropad=InputDict["IsZeropadBeforeStack"])

		C.misc["reference"] = Ctemp.corr[:,1]

	elseif InputDict["IsAlternateRefChannel"]
		#===NOTE===
		Using alternative channel as reference if the identical reference is not found.
		# e.g. When the reference duration contains only BP.LCCB..BP1-BP.MMNB..BP1,
		and compute short-time stack with BP.LCCB..SP1-BP.MMNB..SP1, use BP.LCCB..BP1-BP.MMNB..BP1
		as reference.
		The station keys are reformatted from BP.LCCB..BP1-BP.MMNB..BP1 to BP.LCCB-BP.MMNB-11,
		then search the reference with findfirst().
		=========#

		# refkeypairs = [(original key, no-channel key)]
		refkeypairs = Tuple[]
		for key in keys(ReferenceDict)
			originalstationkey, freqkey = split(key, "/")
			sta1, sta2 = split(originalstationkey, "-")
			net1, sta1, loc1, cha1 = split(sta1, ".")
			net2, sta2, loc2, cha2 = split(sta2, ".")
			nochankey = joinpath("$(net1).$(sta1)-$(net2).$(sta2)-$(cha1[end])$(cha2[end])", freqkey)
			push!(refkeypairs, (key, nochankey))
		end

		# reformat request refdictpath
		originalstationkey, freqkey = split(refdictpath, "/")
		sta1, sta2 = split(originalstationkey, "-")
		net1, sta1, loc1, cha1 = split(sta1, ".")
		net2, sta2, loc2, cha2 = split(sta2, ".")
		nochan_refdictpath = joinpath("$(net1).$(sta1)-$(net2).$(sta2)-$(cha1[end])$(cha2[end])", freqkey)

		alt_refid = findfirst(x -> x== nochan_refdictpath, (y->y[2]).(refkeypairs))
		isnothing(alt_refid) && return nothing # reference is not found
		alt_refdictpath = refkeypairs[alt_refid][1]

		Ctemp = ReferenceDict[alt_refdictpath]

		slice_codawindow!(Ctemp,
							InputDict["background_vel"],
							InputDict["coda_Qinv"],
							InputDict["min_ballistic_twin"],
							InputDict["max_coda_length"],
							attenuation_minthreshold=InputDict["slice_minthreshold"],
							zeropad=InputDict["IsZeropadBeforeStack"])

		C.misc["reference"] = Ctemp.corr[:,1]
	end

    return nothing
end
