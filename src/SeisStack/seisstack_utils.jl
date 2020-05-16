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
	println(keys(ReferenceDict))
	println(refdictpath)

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
	end

    return nothing
end
