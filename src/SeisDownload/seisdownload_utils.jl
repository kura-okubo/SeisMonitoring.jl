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


"""
convert_tmpfile(InputDict::Dict)

convert temporal file in "./seisdownload_tmp" to prescribed format.
It has salvage mode, which allows to compile the temporal files in the case of failing during the download.
"""
function convert_tmpfile(InputDict::OrderedDict; salvage::Bool=false)

    println("-------START CONVERTING--------")
    paths   = SeisIO.ls(InputDict["tmpdir_dl"])
    fodir   = InputDict["fodir"]
    fmt     = InputDict["outputformat"]

    if uppercase(fmt) == "JLD2"
        fopath = joinpath(fodir, "RawData.jld2")
		fo = jldopen(fopath, "w")
    elseif uppercase(fmt) == "ASDF"
        fopath = joinpath(fodir, "RawData.h5")
    else
        @error("outputformat: $(outputformat) is not available (JLD2 or ASDF).")
    end

    varnamelist     = []
	stationlist		= []

    @simd for path in paths
        #println(path)
        S = try
            rseis(path)[1]
        catch y
            #println(y)
        end

        for ii = 1:S.n #loop at each seis channel

            # make station list
            staid = S[ii].id
            if isempty(filter(x -> x==staid, stationlist))
                push!(stationlist, staid)
            end

            # save data (checking whether it's already in the jld2 because it causes an error)
            #parse info
            s_str = string(u2d(S[ii].t[1,2]*1e-6))[1:19]

            #NOTE: if time id is not matched with DLtimestamplist, this download is discarded for consistency even if it has some data.
            #Here we also indirectly check if download_margin is correctly manipulated.
            if !isempty(filter(x -> x==s_str[1:19], InputDict["starttimelist"]))
				# compute end time:
				et = S[ii].t[1,2]*1e-6 + (S[ii].t[end,1]-1)/S[ii].fs
				e_str=string(u2d(et))[1:19]

				# split network, station, location and channel
				net, sta, loc, cha = split(S[ii].id, ".")

				groupname = joinpath("Waveforms", join([net, sta], "."))
				varname	  = join([net, sta, loc, cha], ".")*"__"*s_str*"__"*e_str*"__"*lowercase(cha)

                # select output format

				if isempty(filter(x -> x==varname, varnamelist))
					push!(varnamelist, varname)
                	if fmt == "JLD2"
                        fo[joinpath(groupname,varname)] = S[ii]
                	elseif fmt == "ASDF"
                    	write_hdf5(fopath, S[ii], add=true)
					end
				end
            end
        end

        if !InputDict["Istmpfilepreserved"]
			rm(path)
		end
    end

    if fmt == "JLD2"
        JLD2.close(fo)
	end

    return nothing
end
