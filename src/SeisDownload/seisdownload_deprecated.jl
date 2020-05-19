# # Check maximum memory allocation
# if sizeof(Stemp)/1024/1024/1024 > InputDict["MAX_MEM_PER_CPU"]
# 	@warn "maximam allocation of memory per cpu exceeds predescribed MAX_MEM_PER_CPU.
# 	This may cause transient memory leak, so please track the memory usage." AllocatedMemory_GB=sizeof(Stemp)/1024/1024/1024
# end



# """
#     seisdownload_EARTHQUAKE(startid, InputDict::Dict)
#
# Download seismic data, removing instrumental response and saving into JLD2 file.
#
# # Arguments
# - `startid`         : start time id in starttimelist
# - `InputDict::Dict` : dictionary which contains request information
# """
# function seisdownload_EARTHQUAKE(startid, InputDict::Dict)
#
# 	method		    = InputDict["method"]
# 	event		    = InputDict["event"]
# 	reg			    = InputDict["reg"]
#     pre_filt        = InputDict["pre_filt"]
#
#     #show progress
#     if mod(startid, round(0.1*length(event))+1) == 0
#         println("start downloading event number: $startid")
#     end
#     S = SeisData()
#
#     #---download data---#
# 	for j = 1:length(event[startid]["pickphase"])
# 		net = event[startid]["pickphase"][j]["net"]
# 		sta = event[startid]["pickphase"][j]["sta"]
# 		loc = event[startid]["pickphase"][j]["loc"]
# 		cha = event[startid]["pickphase"][j]["cha"]
# 		src = event[startid]["pickphase"][j]["src"]
#
# 		starttime = string(event[startid]["pickphase"][j]["starttime"])
# 		endtime = string(event[startid]["pickphase"][j]["endtime"])
#
# 		# make multiple request str
#
#     	requeststr = join([net,sta,loc,cha], ".")
#
# 		if InputDict["IsLocationBox"]
# 			# request with lat-lon box
#     		#argv = [method, requeststr, starttime, endtime, InputDict["reg"], 0, src, false, "$requeststr.$startid.xml"]
# 		    ex = :(get_data($(method), $(requeststr), s=$(starttime), t=$(endtime), reg=$(InputDict["reg"]), v=$(0), src=$(src), xf=$("$requeststr.$startid.xml")))
# 		    Stemp = check_and_get_data(ex, requeststr)
# 		else
# 			# request with lat-lon box
#     		#argv = [method, requeststr, starttime, endtime, 0, src, false, "$requeststr.$startid.xml"]
# 			ex = :(get_data($(method), $(requeststr), s=$(starttime), t=$(endtime), v=$(0), src=$(src), xf=$("$requeststr.$startid.xml")))
# 		    Stemp = check_and_get_data(ex, requeststr)
# 		end
#
# 	    if Stemp.misc[1]["dlerror"] == 0 && InputDict["IsResponseRemove"]
# 			if InputDict["IsXMLfileRemoved"]
# 				rm("$requeststr.$startid.xml")
# 			else
# 				mkpath("./stationxml")
# 				mv("$requeststr.$startid.xml", "./stationxml/$requeststr.$startid.xml", force=true)
# 			end
# 		else
# 			if InputDict["IsXMLfileRemoved"]
# 				rm("$requeststr.$startid.xml")
# 			else
# 				mkpath("./stationxml")
# 				mv("$requeststr.$startid.xml", "./stationxml/$requeststr.$startid.xml", force=true)
# 			end
# 	    end
#
# 		#fill gap with zero
# 		SeisIO.ungap!(Stemp, m=true)
# 		replace!(Stemp.x, NaN=>0)
#
# 		append!(S, Stemp)
# 	end
#
#     return S
# end


# method		    = InputDict["method"]
# event		    = InputDict["event"]
# reg			    = InputDict["reg"]
# fopath          = InputDict["fopath"]
#
# #save info into jld2
# jldopen(fopath, "w") do file
# 	file["info/method"]  = method;
# 	file["info/event"]   = event;
# 	file["info/reg"]     = reg
# 	file["info/fopath"]  = fopath
# end
#
# #Test download to evaluate use of memory and estimate download time.
# max_num_of_processes_per_parallelcycle = testdownload(InputDict, length(event), MAX_MEM_PER_CPU)
#
# if max_num_of_processes_per_parallelcycle < 1
# 	error("Memory allocation is not enought (currently $MAX_MEM_PER_CPU [GB]). Please inclease MAX_MEM_PER_CPU or decrease number of stations")
# end
#
# if max_num_of_processes_per_parallelcycle >= length(event)
#
# 	S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), 1:length(event))
#
# 	# save data to jld2
# 	file = jldopen(fopath, "r+")
# 	unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")
#
# 	for ii = 1:size(S)[1] #loop at each starttime
# 		varname = joinpath("event",InputDict["event"][ii]["origin"]["time"][1:end-4])
# 		file[varname] = S[ii]
# 	end
#
# 	JLD2.close(file)
# 	JLD2.close(unavalilablefile)
#
# else
#
# 	#parallelization by time
# 	pitr = 1
#
# 	while pitr <= length(event)
#
# 		startid1 = pitr
# 		startid2 = pitr + max_num_of_processes_per_parallelcycle - 1
#
# 		if startid1 == length(event)
# 			#use one
# 			S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid1)
#
# 		elseif startid2 <= length(event)
# 			#use all processors
# 			S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid2)
#
# 		else
# 			#use part of processors
# 			startid2 = startid1 + mod(length(event), max_num_of_processes_per_parallelcycle) - 1
# 			S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid2)
# 		end
#
# 		# save data to jld2
# 		file = jldopen(fopath, "r+")
# 		unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")
#
#
# 		for ii = 1:size(S)[1] #loop at each starttime
# 			varname = joinpath("event",InputDict["event"][startid1+ii-1]["origin"]["time"][1:end-4])
# 			file[varname] = S[ii]
# 		end
#
# 		JLD2.close(file)
# 		JLD2.close(unavalilablefile)
#
# 		pitr += max_num_of_processes_per_parallelcycle
#
# 		#println("pitr: $pitr")
# 	end
# end

# """
#     get_timestamplist(stlist::Array{DateTime, 1})
#
#     returns list of timestamp: Format = "Year_Julianday_Starttime".
#
# """
# function get_timestamplist(stlist::Array{Any,1})
#
#     timestamplist = []
#     for stid = 1:length(stlist)
#         yj = parse(Int64,stlist[stid][1:4])
#         dj = md2j(yj, parse(Int64,stlist[stid][6:7]), parse(Int64,stlist[stid][9:10]))
#         groupname    = string(yj)*"."*string(dj)*"."*stlist[stid][11:19] #Year_Julianday_Starttime
#         push!(timestamplist, groupname)
#     end
#
#     return timestamplist
#
# end

# """
#     get_stationlist(network::Array{String, 1}, station::Array{String, 1}, location::Array{String, 1}, channel::Array{String, 1})
#
#     returns list of request strings:
#
# """
# function get_stationlist(network::Array{String, 1}, station::Array{String, 1}, location::Array{String, 1}, channel::Array{String, 1})
#
#     stationlist = []
#     for networkid = 1:length(network)
#         for stationid = 1:length(station)
#             for locationid = 1:length(location)
#                 for channelid = 1:length(channel)
#                     requeststr = @sprintf("%s.%s.%s.%s", network[networkid], station[stationid], location[locationid], channel[channelid])
#                     push!(stationlist, requeststr)
#                 end
#             end
#         end
#     end
#
#     return stationlist
#
# end

# """
#     testdownload(InputDict::Dict{String,Any} numofitr::Int64)
#
#     print stats of download and return max_num_of_processes_per_parallelcycle
#
# # Output
#  -`max_num_of_processes_per_parallelcycle`: maximum number of processes for one request
#
# Deprecated
# """
# function testdownload(InputDict::Dict{String,Any}, numofitr::Int64)
#
#     DownloadType    = InputDict["DownloadType"]
#
#     trial_id          = 1
#     InputDict_test = deepcopy(InputDict)
#     test_suceededflag = false
#     println("-------TEST DOWNLOAD START-----------")
#
#     if DownloadType == "Noise" || DownloadType == "noise"
#
#         while !test_suceededflag && trial_id < length(InputDict["starttimelist"])
#             # select test request
#             for j = 1:length(InputDict["stationinfo"]["stationlist"])
#                 InputDict_test["stationinfo"]["stationlist"] = [InputDict["stationinfo"]["stationlist"][j]]
#
#                 global t1 = @elapsed global dlerror = seisdownload_NOISE(trial_id, InputDict_test, testdownload=true) #[s]
#
#                 dl = [dlerror[i] for i in 1:length(dlerror)]
#                 if issubset(0, dl)
#                     test_suceededflag = true
#                     break;
#                 end
#             end
#             trial_id += 1
#         end
#
#     elseif  DownloadType == "Earthquake" || DownloadType == "earthquake"
#
#         while !test_suceededflag && trial_id < length(InputDict["starttimelist"])
#             # select test request
#             for j = 1:length(InputDict["stationinfo"]["stationlist"])
#                 InputDict_test["stationinfo"]["stationlist"] = InputDict["stationinfo"]["stationlist"][j]
#
#                 global t1 = @elapsed global dlerror = seisdownload_EARTHQUAKE(trial_id, InputDict_test) #[s]
#
#                 dl = [dlerror[i] for i in 1:length(dlerror)]
#                 if issubset(0, dl)
#                     test_suceededflag = true
#                     break;
#                 end
#             end
#             trial_id += 1
#         end
#     end
#
#     if !test_suceededflag
#         error("All requests you submitted with input dictionary was failed. Please check the station availability in your request.")
#     end
#
#     estimated_downloadtime = now() + Second(round(3 * t1 * length(InputDict["DLtimestamplist"]) * numofitr / nprocs()))
#
#     #println(mem_per_requestid)
#     #println(max_num_of_processes_per_parallelcycle)
#     println("-------DOWNLOAD STATS SUMMARY--------")
#
#     println(@sprintf("Number of processes is %d.", nprocs()))
#
#     # evaluate total download size by searching tmp directory
#     tmppath = InputDict["tmppath"]
#     s = read(`du -s -k $tmppath`, String)
#     hdduse = parse(Int, split(s)[1])
#
#     totaldownloadsize = hdduse * numofitr * length(InputDict["stationinfo"]["stationlist"])
#     if totaldownloadsize < 1024 * 1024 # less than 1 GB
#         totaldownloadsize = totaldownloadsize / 1024 #[MB]
#         sizeunit = "MB"
#     else
#         totaldownloadsize = totaldownloadsize / 1024 / 1024
#         sizeunit = "GB"
#     end
#
#     println(@sprintf("Total download size will be %4.2f [%s].", 0.8 * totaldownloadsize, sizeunit)) #0.8: considering compression efficiency
#     println(@sprintf("Download will finish at %s.", round(estimated_downloadtime, Dates.Second(1))))
#     println("*We have a time lag with downloading time above, like in 10 minutes or so.*")
#     println("*This estimation also changes if some download requests fail and are skipped.*")
#     println("-------START DOWNLOADING-------------")
#
#     return nothing
#
# end
#
# """
# printparams(param::Dict)
#
# print parameters
# """
# function printparams(param::Dict)
#     printstyled("-----------Input Parameters-----------\n"; color=:cyan, bold=true)
#     for key in keys(param)
#         if length(string(param["$key"])) > 60
#             param_str = string(param["$key"])[1:30]*"..."
#         else
#             param_str = string(param["$key"])
#         end
#         println(@sprintf("%-24s = %-10s", key, param_str))
#     end
# end


# """
# convert_tmpfile(InputDict::Dict)
#
# convert temporal file in "./seisdownload_tmp" to prescribed format.
# It has salvage mode, which allows to compile the temporal files in the case of failing during the download.
# """
# function convert_tmpfile(InputDict::OrderedDict; salvage::Bool=false)
#
#     paths   = SeisIO.ls(InputDict["tmpdir_dl"])
#     fodir   = InputDict["fodir"]
#     fmt     = InputDict["outputformat"]
#
#     if uppercase(fmt) == "JLD2"
#         fopath = joinpath(fodir, "RawData.jld2")
# 		fo = jldopen(fopath, "w")
#     elseif uppercase(fmt) == "ASDF"
#         fopath = joinpath(fodir, "RawData.h5")
#     else
#         @error("outputformat: $(outputformat) is not available (JLD2 or ASDF).")
#     end
#
#     varnamelist     = []
# 	stationlist		= []
#
#     for path in paths
#         #println(path)
#         S = try
#             rseis(path)[1]
#         catch y
#             #println(y)
#         end
#
#         for ii = 1:S.n #loop at each seis channel
#
#             # make station list
#             staid = S[ii].id
#             if isempty(filter(x -> x==staid, stationlist))
#                 push!(stationlist, staid)
#             end
#
#             # save data (checking whether it's already in the jld2 because it causes an error)
#             #parse info
# 			isempty(S[ii].t) && continue;
#             s_str = string(u2d(S[ii].t[1,2]*1e-6))[1:19]
#
#             #NOTE: if time id is not matched with DLtimestamplist, this download is discarded for consistency even if it has some data.
#             #Here we also indirectly check if download_margin is correctly manipulated.
#             if !isempty(filter(x -> x==s_str[1:19], InputDict["starttimelist"]))
# 				# compute end time:
# 				et = S[ii].t[1,2]*1e-6 + (S[ii].t[end,1]-1)/S[ii].fs
# 				e_str=string(u2d(et))[1:19]
#
# 				# split network, station, location and channel
# 				net, sta, loc, cha = split(S[ii].id, ".")
#
# 				groupname = joinpath("Waveforms", join([net, sta], "."))
# 				varname	  = join([net, sta, loc, cha], ".")*"__"*s_str*"__"*e_str*"__"*lowercase(cha)
#
#                 # select output format
#
# 				if isempty(filter(x -> x==varname, varnamelist))
# 					push!(varnamelist, varname)
#                 	if fmt == "JLD2"
#                         fo[joinpath(groupname,varname)] = S[ii]
#                 	elseif fmt == "ASDF"
#                     	write_hdf5(fopath, S[ii], add=true)
# 					end
# 				end
#             end
#         end
#
#         if !InputDict["Istmpfilepreserved"]
# 			rm(path)
# 		end
#     end
#
#     if fmt == "JLD2"
#         JLD2.close(fo)
# 	end
#
#     return nothing
# end
