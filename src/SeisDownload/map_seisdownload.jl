"""
    map_seisdownload_NOISE(startid, InputDict::Dict)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `startid`         : start time id in starttimelist
- `InputDict::Dict` : dictionary which contains request information
"""
function map_seisdownload_NOISE(startid, InputDict::OrderedDict; testdownload::Bool=false)

	StationDataFrame = jldopen(InputDict["requeststation_file"]) # (e.g. requeststation.jld2 contains Dataframe)

    download_time_unit = InputDict["download_time_unit"]

	fodir 				 = InputDict["fodir"]
	tmpdir 				 = InputDict["tmpdir_dl"]
	stationxml_dir		 = InputDict["stationxml_dir"]
	requeststation_file	 = 	InputDict["requeststation_file"]

	#SeisIO getdata option
	if !haskey(InputDict, "get_data_opt")
		InputDict["get_data_opt"] = [true, true, true, true, true]#: [unscale, demean, detrend, taper, ungap]
	end

	if !haskey(InputDict, "sampling_frequency")
		InputDict["sampling_frequency"] = false # default margin: 5 minutes
	end

	if !haskey(InputDict, "download_margin")
		InputDict["download_margin"] = 5 * 60 # default margin: 5 minutes
	end

    #make stlist at all processors
    starttimelist = InputDict["starttimelist"]

    #show progress
    if starttimelist[startid][end-8:end] == "T00:00:00" && !testdownload
        println("start downloading $(starttimelist[startid])")
	end

	# import request_networkchanks (see )
	if !ispath(requeststation_file)
		error("request station file: $(requeststation_file) is not found.")
	end

	request_src_chanks = jldopen(requeststation_file, "r")

    for src in keys(request_src_chanks)
        #---download data---#

		requeststrs = get_requeststr(request_src_chanks[src], InputDict["numstationperrequest"])

		for requeststr in requeststrs
			#NOTE: method is currently fixed "FDSN", which downloads station xml and waveforms.
			method = "FDSN"

			# including download margin
			starttime = string(DateTime(starttimelist[startid]) - Second(InputDict["download_margin"]))
			dltime = download_time_unit + 2 * InputDict["download_margin"]

			stationxml_path = joinpath(stationxml_dir, "$src.$starttime.xml")

			if InputDict["IsLocationBox"]
		        ex = :(get_data($(method), $(requeststr), s=$(starttime), t=$(dltime), reg=$(InputDict["reg"]),
				 v=$(0), src=$(src), xf=$(stationxml_path), unscale=$(InputDict["get_data_opt"][1]),
				  demean=$(InputDict["get_data_opt"][2]), detrend=$(InputDict["get_data_opt"][3]),taper=$(InputDict["get_data_opt"][4]),
				  ungap=$(InputDict["get_data_opt"][5]), rr=$(InputDict["IsResponseRemove"])))

		        t_dl = @elapsed Stemp = check_and_get_data(ex, requeststr)
			else
				ex = :(get_data($(method), $(requeststr), s=$(starttime), t=$(dltime),
				 v=$(0), src=$(src), xf=$(stationxml_path),unscale=$(InputDict["get_data_opt"][1]),
				  demean=$(InputDict["get_data_opt"][2]), detrend=$(InputDict["get_data_opt"][3]),taper=$(InputDict["get_data_opt"][4]),
				  ungap=$(InputDict["get_data_opt"][5]), rr=$(InputDict["IsResponseRemove"])))

				t_dl = @elapsed Stemp = check_and_get_data(ex, requeststr)
			end

			if Stemp == 1
				# error this entire request
				return 1
			end

			Isdataflag = false
			# manipulate download_margin
			# manipulate_tmatrix!(Stemp, starttime, InputDict)
			# NOTE: using SeisIO.sync()
			stsync = DateTime(starttimelist[startid])
			etsync = u2d(d2u(stsync) + download_time_unit)
			# println([stsync, etsync])
			SeisIO.sync!(Stemp, s=stsync, t=etsync, v=0)

			if InputDict["sampling_frequency"] isa Number
				for j = 1:Stemp.n
					if Stemp.misc[j]["dlerror"] == 0 && !isempty(Stemp[j].t)
						# downsample
						if Stemp.fs[j] > InputDict["sampling_frequency"]
							SeisIO.resample!(Stemp, chans=j, fs=float(InputDict["sampling_frequency"]))
						end
					end
				end
			end

			# Save temp file
			requeststr_all = join(requeststr, "--")
			fname_out = join([requeststr_all, string(stsync), string(etsync), src], "__")*".dat"
			# fname_out = join([String(y),
			# 			string(j),
			# 			replace(split(starttimelist[startid], 'T')[2], ':' => '.'),
			# 			requeststr,
			# 			"FDSNWS",
			# 			src,"dat"],
			# 			'.')
			# save as intermediate binary file
			t_write = @elapsed wseis(joinpath(InputDict["tmpdir_dl"], fname_out), Stemp)

			if !InputDict["IsXMLfilepreserved"] && ispath(stationxml_path)
				rm(stationxml_path)
			end
		end
    end

    return 0

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
function check_and_get_data(ex::Expr, requeststr::Array{String,1})

	# we try the same download request upto 3 times because it sometimes fails
	# due to network error.
	download_itr = 1

	S = SeisData()

	while download_itr <= 3
		S = try
			#remove comment out below if you want to print contents of get_data()
			# println(ex)
			eval(ex);
		catch
		end

		if !isnothing(S) && !isempty(S)
			break;
		else
			download_itr += 1
			#println("retry downloading $(download_itr).")
		end
	end

	if isnothing(S) || isempty(S)
		return 1
	end

	for j = 1:S.n
		!isnothing(S[j]) && !isempty(S[j]) && !isempty(S[j].t) ? S.misc[j]["dlerror"] = 0 : S.misc[j]["dlerror"] = 1
	end

	return S
end

#
# """
#     manipulate_tmatrix!(S::SeisData, InputDict::Dict{String,Any})
#
# manipulate time matrix to remove download margin
# """
# function manipulate_tmatrix!(S::SeisData, starttime::String, InputDict::OrderedDict)
#
#     for i = 1:S.n
#
# 		if S.misc[i]["dlerror"] == 1
# 			continue;
# 		end
#
#         download_margin = InputDict["download_margin"]
#         download_time_unit    = InputDict["download_time_unit"]
#         requeststr = S.id[i]
#
# 		# NOTE: 2020/2/8
# 		# We decided to use data ungap to avoid data inconsitency
# 		# due to data gap and time shift.
# 		if size(S.t[i], 1) > 2
# 			warning("ungap! is not applied yet to SeisData, which may cause
# 			data inconsistency. So applying ungap.")
# 			ungap!(S[i])
# 		end
#
#         tvec = collect(0:S.t[i][end,1]-1) ./ S.fs[i]
#         tlen = trunc(Int, download_time_unit * S.fs[i])
#
#         si = findfirst(x -> tvec[x] >= download_margin, 1:length(tvec))
#
# 		#println([si, download_margin * S.fs[i]])
# 		#println([string(u2d(S.t[i][1,2] * 1e-6))[1:19],starttime])
# 		# Check if data is within request time window AND start time is equal
# 		# at the order of second to what is requested
# 		# rounding subsecond error in downloading
# 		# NOTE: thie is just used to check the consistency between requested
# 		# time and truncated time;
# 		# SeisNoise.phase_shift! is applied when seisxcorr is performed.
#
# 		tsync = round(Int, S.t[i][1,2] * 1e-6) * 1e6
#
# 		if isnothing(si)
# 			#downloaded data is only within margin; nodata in the requested time window
# 			S.misc[i]["dlerror"] = 1
#             S.x[i] = zeros(0)
#
#         elseif si < download_margin * S.fs[i] || string(u2d(tsync * 1e-6))[1:19] != starttime
#             #println("data missing or starttime not match.")
#             S.misc[i]["dlerror"] = 1
#             S.x[i] = zeros(0)
#         else
#             #println("manipulate")
# 			#===
# 			Shift time index of si:
# 			|--------|----------------------------------------|--------|
# 			  margin              requested data (tlen)         margin
# 			        [si]                                 [si]+tlen-1
#
# 			===#
#             ei = trunc(Int, min(si+tlen-1, S.t[i][end,1]))
#             x_shifted = zeros(tlen)
#             copyto!(x_shifted, S.x[i][si:ei])
#
# 			# # manipulate time matrix
# 			# t_shifted = deepcopy(hcat(S.t[i][:,1] .- (si-1), S.t[i][:,2]))
# 			# # remove gap outside of time window
# 			# ri = findall(x -> t_shifted[x, 1] < 1 || t_shifted[x, 1] >  tlen, 1:size(t_shifted, 1))
# 			# t_shifted = t_shifted[setdiff(1:end, ri), :]
# 			t_init = [1 trunc(Int, S.t[i][1,2] + float(download_margin)*1e6)]
# 			t_last = [tlen 0]
# 			t_shifted = vcat(t_init, t_last)
# 			S.x[i] = x_shifted
# 			S.t[i] = t_shifted
#
#         end
#     end
#     #print("after")
#     #println(S)
# end
