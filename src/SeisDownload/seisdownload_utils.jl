using SeisIO, JLD2, Dates
"""

get_starttimelist(st::DateTime, et::DateTime, unittime::Float64)
calculate start time list for parallel downloading

    st: start time
    et: end time
    unittime: unit time in Second

    this function returns
    stlist: list of start time

    e.g.
    st = DateTime(2019,1,1,0,0,0)
    et = DateTime(2019,1,1,12,0,0)
    unittime = 3600

    stlist = get_starttimelist(st, et, unittime)

"""
function get_starttimelist(st::DateTime, et::DateTime, unittime::Real)

    reftime = st
    stlist = []

    while reftime < et
        push!(stlist, string(reftime))
        reftime += Dates.Second(float(unittime))
    end

    return stlist
end


"""
	get_requeststr(df::DataFrame, numstationperrequest::Int))

return request str following web_chanspec of SeisIO.get_data.

- 'numstationperrequest::Int': number of stations per one request to avoid too many request statinos at one time.
"""
function get_requeststr(df::DataFrame, numstationperrequest::Int)

	reqstrs = Array{Array{String,1},1}(undef, 0) # array of request stations per HTTP request
	for inds = Iterators.partition(1:size(df)[1], numstationperrequest)
		reqstr  = String[]
		for i in inds
			rst = join([df.network[i], df.station[i], df.location[i], df.channel[i]], ".")
			push!(reqstr, rst)
		end
		push!(reqstrs, reqstr)
	end
	#
	# for i = 1:size(df)[1]
	# 	rst = join([df.network[i], df.station[i], df.location[i], df.channel[i]], ".")
	# end
	return reqstrs
end
