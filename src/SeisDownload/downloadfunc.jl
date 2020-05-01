module DownloadFunc

using SeisIO, Dates

export seisdownload_NOISE, seisdownload_EARTHQUAKE

"""
    seisdownload_NOISE(startid, InputDict::Dict)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `startid`         : start time id in starttimelist
- `InputDict::Dict` : dictionary which contains request information
"""
function seisdownload_NOISE(startid, InputDict::Dict; testdownload::Bool=false)

    #stationlist
    stationlist     = InputDict["stationinfo"]["stationlist"]
    method      	= InputDict["stationinfo"]["stationmethod"]
    src             = InputDict["stationinfo"]["stationsrc"]
    starttime       = InputDict["starttime"]
    endtime         = InputDict["endtime"]
    DL_time_unit    = InputDict["DL_time_unit"]

	#SeisIO getdata option
	if !haskey(InputDict, "get_data_opt")
		InputDict["get_data_opt"] = [true, true, true, true, true]#: [unscale, demean, detrend, taper, ungap]
	end

	if !haskey(InputDict, "savesamplefreq")
		InputDict["savesamplefreq"] = false # default margin: 5 minutes
	end

	if !haskey(InputDict, "download_margin")
		InputDict["download_margin"] = 5 * 60 # default margin: 5 minutes
	end

    #make stlist at all processors

    starttimelist = InputDict["starttimelist"]

    #show progress
    if starttimelist[startid][end-8:end] == "T00:00:00" && !testdownload
        println("start downloading $(starttimelist[startid])")
	elseif testdownload
		print(".")
    end

	dlerror = []

    for i = 1:length(stationlist)
        #---download data---#
        requeststr = stationlist[i]

		#println(starttimelist[startid])

		# including download margin
		starttime = string(DateTime(starttimelist[startid]) - Second(InputDict["download_margin"]))
		dltime = DL_time_unit + 2 * InputDict["download_margin"]

		#println(starttime)
		#println(dltime)

		if InputDict["IsLocationBox"]

	        ex = :(get_data($(method[i]), $(requeststr), s=$(starttime), t=$(dltime), reg=$(InputDict["reg"]),
			 v=$(0), src=$(src[i]), xf=$("$requeststr.$startid.xml"), unscale=$(InputDict["get_data_opt"][1]),
			  demean=$(InputDict["get_data_opt"][2]), detrend=$(InputDict["get_data_opt"][3]),taper=$(InputDict["get_data_opt"][4]),
			  ungap=$(InputDict["get_data_opt"][5]), rr=$(InputDict["IsResponseRemove"])))

	        t_dl = @elapsed Stemp = check_and_get_data(ex, requeststr)
		else

			ex = :(get_data($(method[i]), $(requeststr), s=$(starttime), t=$(dltime),
			 v=$(0), src=$(src[i]), xf=$("$requeststr.$startid.xml"),unscale=$(InputDict["get_data_opt"][1]),
			  demean=$(InputDict["get_data_opt"][2]), detrend=$(InputDict["get_data_opt"][3]),taper=$(InputDict["get_data_opt"][4]),
			  ungap=$(InputDict["get_data_opt"][5]), rr=$(InputDict["IsResponseRemove"])))

			t_dl = @elapsed Stemp = check_and_get_data(ex, requeststr)
		end

		# Check maximum memory allocation
		if sizeof(Stemp)/1024/1024/1024 > InputDict["MAX_MEM_PER_CPU"]
			@warn "maximam allocation of memory per cpu exceeds predescribed MAX_MEM_PER_CPU.
			This may cause transient memory leak, so please track the memory usage." AllocatedMemory_GB=sizeof(Stemp)/1024/1024/1024
		end

		Isdataflag = false
		# manipulate download_margin
		manipulate_tmatrix!(Stemp, starttime, InputDict)

		for j = 1:Stemp.n
			if Stemp.misc[j]["dlerror"] == 0
				Isdataflag = true
				# downsample
				if InputDict["savesamplefreq"] isa Number
					if Stemp.fs[j] > InputDict["savesamplefreq"]
						SeisIO.resample!(Stemp, chans=j, fs=float(InputDict["savesamplefreq"]))
					end
				end
			end
		end

		# if some of SeisChannels in Stemp have a data, save temp file
		if Isdataflag
			ymd = split(starttimelist[startid], r"[A-Z]")
			(y, m, d) = split(ymd[1], "-")
			j = md2j(y, m, d)
			fname_out = join([String(y),
						string(j),
						replace(split(starttimelist[startid], 'T')[2], ':' => '.'),requeststr,
						"FDSNWS",
						src[i],"dat"],
						'.')

			# save as intermediate binary file
			t_write = @elapsed wseis(InputDict["tmppath"]*"/"*fname_out, Stemp)
		end


		if InputDict["IsXMLfileRemoved"] && ispath("$requeststr.$startid.xml")
			rm("$requeststr.$startid.xml")
		end

		push!(dlerror, !Isdataflag)

		#print("[dltime, wtime, fraction of writing]: ")
		#println([t_dl, t_write, t_write/(t_dl+t_write)])
    end

    return dlerror
end



"""
    seisdownload_EARTHQUAKE(startid, InputDict::Dict)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `startid`         : start time id in starttimelist
- `InputDict::Dict` : dictionary which contains request information
"""
function seisdownload_EARTHQUAKE(startid, InputDict::Dict)

	method		    = InputDict["method"]
	event		    = InputDict["event"]
	reg			    = InputDict["reg"]
    pre_filt        = InputDict["pre_filt"]

    #show progress
    if mod(startid, round(0.1*length(event))+1) == 0
        println("start downloading event number: $startid")
    end
    S = SeisData()

    #---download data---#
	for j = 1:length(event[startid]["pickphase"])
		net = event[startid]["pickphase"][j]["net"]
		sta = event[startid]["pickphase"][j]["sta"]
		loc = event[startid]["pickphase"][j]["loc"]
		cha = event[startid]["pickphase"][j]["cha"]
		src = event[startid]["pickphase"][j]["src"]

		starttime = string(event[startid]["pickphase"][j]["starttime"])
		endtime = string(event[startid]["pickphase"][j]["endtime"])

		# make multiple request str

    	requeststr = join([net,sta,loc,cha], ".")

		if InputDict["IsLocationBox"]
			# request with lat-lon box
    		#argv = [method, requeststr, starttime, endtime, InputDict["reg"], 0, src, false, "$requeststr.$startid.xml"]
		    ex = :(get_data($(method), $(requeststr), s=$(starttime), t=$(endtime), reg=$(InputDict["reg"]), v=$(0), src=$(src), xf=$("$requeststr.$startid.xml")))
		    Stemp = check_and_get_data(ex, requeststr)
		else
			# request with lat-lon box
    		#argv = [method, requeststr, starttime, endtime, 0, src, false, "$requeststr.$startid.xml"]
			ex = :(get_data($(method), $(requeststr), s=$(starttime), t=$(endtime), v=$(0), src=$(src), xf=$("$requeststr.$startid.xml")))
		    Stemp = check_and_get_data(ex, requeststr)
		end

	    if Stemp.misc[1]["dlerror"] == 0 && InputDict["IsResponseRemove"]
	        #Remove_response_obspy.remove_response_obspy!(Stemp, "$requeststr.$startid.xml", pre_filt=pre_filt, zeropadlen = float(30*60), output="VEL")
			if InputDict["IsXMLfileRemoved"]
				rm("$requeststr.$startid.xml")
			else
				mkpath("./stationxml")
				mv("$requeststr.$startid.xml", "./stationxml/$requeststr.$startid.xml", force=true)
			end
		else
			if InputDict["IsXMLfileRemoved"]
				rm("$requeststr.$startid.xml")
			else
				mkpath("./stationxml")
				mv("$requeststr.$startid.xml", "./stationxml/$requeststr.$startid.xml", force=true)
			end
	    end

		#fill gap with zero
		SeisIO.ungap!(Stemp, m=true)
		replace!(Stemp.x, NaN=>0)

		append!(S, Stemp)
	end

    return S
end


"""
    check_and_get_data(ex::Expr, requeststr::String)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `ex::Expr`        : expression of get data includin all request information

# Output
- `S::SeisData`     : downloaded SeisData
- `requeststr::String`     : request channel (e.g. "BP.LCCB..BP1")
"""
function check_and_get_data(ex::Expr, requeststr::String)
	try
		#comment out below if you want to print contents of get_data()
		#println(ex)
		S = eval(ex);

		for j = 1:S.n
			if !isempty(S.x[j])
				#download succeeded
				S.misc[j]["dlerror"] = 0
			else
				S.misc[j]["dlerror"] = 1
			end
		end

		return S

	catch y
		#println(y)
		S = SeisData(1)
		S.misc[1]["dlerror"] = 1
		S.id[1] = requeststr
		note!(S, 1, "station is not available for this request.")
		return S
	end
end


"""
    manipulate_tmatrix!(S::SeisData, InputDict::Dict{String,Any})

    manipulate time matrix to remove download margin

"""

function manipulate_tmatrix!(S::SeisData, starttime::String, InputDict::Dict{String,Any})

    #print("before")
    #println(S)
    for i = 1:S.n

		if S.misc[i]["dlerror"] == 1
			continue;
		end

        download_margin = InputDict["download_margin"]
        DL_time_unit    = InputDict["DL_time_unit"]
        #requeststr = "NC.PLO..EHZ" #NC.PCC..EHZ
        requeststr = S.id[i] #NC.PCC..EHZ

		# NOTE: 2020/2/8
		# We decided to use data ungap to avoid data inconsitency
		# due to data gap and time shift.
		if size(S.t[i], 1) > 2
			warning("ungap! is not applied yet to SeisData, which may cause
			data inconsistency. So applying ungap.")
			ungap!(S[i])
		end

        tvec = collect(0:S.t[i][end,1]-1) ./ S.fs[i]
        tlen = trunc(Int, DL_time_unit * S.fs[i])

        si = findfirst(x -> tvec[x] >= download_margin, 1:length(tvec))

		#println([si, download_margin * S.fs[i]])
		#println([string(u2d(S.t[i][1,2] * 1e-6))[1:19],starttime])
		# check if data is within request time window AND start time is equal
		# at the order of second to what is requested
		# rounding subsecond error in downloading
		# NOTE: thie is just used to check the consistency between requested
		# time and truncated time;
		# SeisNoise.phase_shift! is applied when seisxcorr is performed.

		tsync = round(Int, S.t[i][1,2] * 1e-6) * 1e6

		if isnothing(si)
			#downloaded data is only within margin; nodata in the requested time window
			S.misc[i]["dlerror"] = 1
            S.x[i] = zeros(0)

        elseif si < download_margin * S.fs[i] || string(u2d(tsync * 1e-6))[1:19] != starttime
            #println("data missing or starttime not match.")
            S.misc[i]["dlerror"] = 1
            S.x[i] = zeros(0)
        else
            #println("manipulate")
			#===
			Shift time index of si:
			|--------|----------------------------------------|--------|
			  margin              requested data (tlen)         margin
			        [si]                                 [si]+tlen-1

			===#
            ei = trunc(Int, min(si+tlen-1, S.t[i][end,1]))
            x_shifted = zeros(tlen)
            copyto!(x_shifted, S.x[i][si:ei])

			# # manipulate time matrix
			# t_shifted = deepcopy(hcat(S.t[i][:,1] .- (si-1), S.t[i][:,2]))
			# # remove gap outside of time window
			# ri = findall(x -> t_shifted[x, 1] < 1 || t_shifted[x, 1] >  tlen, 1:size(t_shifted, 1))
			# t_shifted = t_shifted[setdiff(1:end, ri), :]
			t_init = [1 trunc(Int, S.t[i][1,2] + float(download_margin)*1e6)]
			t_last = [tlen 0]
			t_shifted = vcat(t_init, t_last)
			S.x[i] = x_shifted
			S.t[i] = t_shifted

        end
    end
    #print("after")
    #println(S)
end


end
