include("utils.jl")
include("downloadfunc.jl")

using .Utils
using .DownloadFunc

using SeisIO, Dates, Printf, JLD2, FileIO, Distributed

"""
    seisdownload(InputDict::Dict)

    Request seismic data and save into jld2.
# Arguments
- `InputDict`    : dictionary which contains request information
"""
function seisdownload(InputDict::Dict)

    Utils.initlogo()

	printparams(InputDict)

	DownloadType    = InputDict["DownloadType"]

	fodir = ""
	sp = splitpath(InputDict["fopath"])
	for i = 1:length(sp)-1
		fodir = joinpath(fodir, sp[i])
	end
	tmppath = joinpath(fodir, "./seisdownload_tmp")
	InputDict["tmppath"] = tmppath
	mkpath(tmppath)

    if DownloadType == "Noise" || DownloadType == "noise"

		#stationlist
		stationlist     = InputDict["stationinfo"]["stationlist"]
		starttime       = InputDict["starttime"]
		endtime         = InputDict["endtime"]
		DL_time_unit    = InputDict["DL_time_unit"]
		DownloadType    = InputDict["DownloadType"]
		fopath          = InputDict["fopath"]

		if mod((endtime - starttime).value,  DL_time_unit) != 0 || (endtime - starttime).value < DL_time_unit
			error("Total download time cannot be devided by Download Time unit; this may cause unexpected result. Abort.")
		end

		# calculate start time list (starttimelist) with each Donwload_time_unit
		starttimelist = Utils.get_starttimelist(starttime, endtime, DL_time_unit)
		# generate DLtimestamplist and ststationlist
		DLtimestamplist = Utils.get_timestamplist(starttimelist)

		InputDict["starttimelist"] = starttimelist
		InputDict["DLtimestamplist"] = DLtimestamplist

		#----Restrict number of processors------#
		#NEVER CHANGE THIS THRESHOLD OTHERWISE IT OVERLOADS THE DATA SERVER
		np = nprocs()
		if np > 100 throw(DomainError(np, "np must be smaller than 100.")) end
		#---------------------------------------#

        # Test download to evaluate use of memory and estimate download time.
		InputDict_test = deepcopy(InputDict) # to avoid overwriting InputDict; unknown bug while deepcopying in testdownload function
		testdownload(InputDict_test, length(starttimelist))

		# Start downloading data
		t_download = @elapsed pmap(x -> seisdownload_NOISE(x, InputDict), 1:length(starttimelist))

		# convert intermediate file to prescibed file format (JLD2, ASDF, ...)
		t_convert = @elapsed convert_tmpfile(InputDict)

		println("---Summary of computational time---")
		println(@sprintf("Total download time:%8.4f[s]", t_download))
		println(@sprintf("Total convert time:%8.4f[s]", t_convert))

		if !InputDict["Istmpfilepreserved"]
			rm(tmppath, recursive=true, force=true)
		end

    elseif  DownloadType == "Earthquake" || DownloadType == "earthquake"

		method		    = InputDict["method"]
		event		    = InputDict["event"]
		reg			    = InputDict["reg"]
		fopath          = InputDict["fopath"]

		#save info into jld2
		jldopen(fopath, "w") do file
			file["info/method"]  = method;
			file["info/event"]   = event;
			file["info/reg"]     = reg
			file["info/fopath"]  = fopath
		end

		#Test download to evaluate use of memory and estimate download time.
		max_num_of_processes_per_parallelcycle = testdownload(InputDict, length(event), MAX_MEM_PER_CPU)

		if max_num_of_processes_per_parallelcycle < 1
			error("Memory allocation is not enought (currently $MAX_MEM_PER_CPU [GB]). Please inclease MAX_MEM_PER_CPU or decrease number of stations")
		end

		if max_num_of_processes_per_parallelcycle >= length(event)

			S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), 1:length(event))

			# save data to jld2
			file = jldopen(fopath, "r+")
			unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")

			for ii = 1:size(S)[1] #loop at each starttime
				varname = joinpath("event",InputDict["event"][ii]["origin"]["time"][1:end-4])
				file[varname] = S[ii]
			end

			JLD2.close(file)
			JLD2.close(unavalilablefile)

		else

			#parallelization by time
			pitr = 1

			while pitr <= length(event)

				startid1 = pitr
				startid2 = pitr + max_num_of_processes_per_parallelcycle - 1

				if startid1 == length(event)
					#use one
					S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid1)

				elseif startid2 <= length(event)
					#use all processors
					S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid2)

				else
					#use part of processors
					startid2 = startid1 + mod(length(event), max_num_of_processes_per_parallelcycle) - 1
					S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid2)
				end

				# save data to jld2
				file = jldopen(fopath, "r+")
				unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")


				for ii = 1:size(S)[1] #loop at each starttime
					varname = joinpath("event",InputDict["event"][startid1+ii-1]["origin"]["time"][1:end-4])
					file[varname] = S[ii]
				end

				JLD2.close(file)
				JLD2.close(unavalilablefile)

				pitr += max_num_of_processes_per_parallelcycle

				#println("pitr: $pitr")
			end
		end


    else
        println("Download type is not known (chose Noise or Earthquake).")
    end

    println("Downloading and Saving data is successfully done.\njob ended at "*string(now()))
    return 0

end
