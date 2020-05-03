include("utils.jl")
include("downloadfunc.jl")

using SeisIO, Dates, Printf, JLD2, FileIO, Distributed

"""
    seisdownload(InputDict::Dict)

    Request seismic data and save into jld2.
# Arguments
- `InputDict`    : dictionary which contains request information
"""
function seisdownload(InputDict::Dict)

	downloadtype 			= InputDict["downloadtype"]
	project_output_dir		= abspath(InputDict["project_output_dir"])
	InputDict["fodir"] 		= joinpath(project_output_dir, "seismicdata")
	tmpdir 					= joinpath(project_output_dir, "seismicdata", "seisdownload_tmp")
	InputDict["tmpdir_dl"] 	= tmpdir
	mkdir(tmpdir)

    if downloadtype == "Noise" || downloadtype == "noise"

		#stationlist
		# stationlist     = InputDict["stationinfo"]["stationlist"]
		starttime       	= InputDict["starttime"]
		endtime         	= InputDict["endtime"]
		download_time_unit  = InputDict["download_time_unit"]

		# Check if total download time can be devided by download time unit.
		if mod((endtime - starttime).value,  download_time_unit) != 0 || (endtime - starttime).value < download_time_unit
			@error("Total download time cannot be devided by Download Time unit; this may cause unexpected result.")
			return 1
		end

		# calculate start time list (starttimelist) with each Donwload_time_unit
		starttimelist = get_starttimelist(starttime, endtime, download_time_unit)
		# generate DLtimestamplist and ststationlist
		# DLtimestamplist = Utils.get_timestamplist(starttimelist)

		InputDict["starttimelist"] = starttimelist
		# InputDict["DLtimestamplist"] = DLtimestamplist

		#----Restrict number of processors------#
		#NEVER CHANGE THIS THRESHOLD OTHERWISE IT OVERLOADS THE DATA SERVER
		np = nprocs()
		if np > 100 throw(DomainError(np, "np must be smaller than 100.")) end
		#---------------------------------------#

        # # Test download to evaluate use of memory and estimate download time.
		# InputDict_test = deepcopy(InputDict) # to avoid overwriting InputDict; unknown bug while deepcopying in testdownload function
		# testdownload(InputDict_test, length(starttimelist))

		# Start downloading data
		println("-------START Downloading--------")

		t_download = @elapsed pmap(x -> seisdownload_NOISE(x, InputDict), 1:length(starttimelist))

		# convert intermediate file to prescibed file format (JLD2, ASDF, ...)
		t_convert = @elapsed convert_tmpfile(InputDict)

		println("---Summary of computational time---")
		println(@sprintf("Total download time:%8.4f[s]", t_download))
		println(@sprintf("Total convert time:%8.4f[s]", t_convert))

		if !InputDict["Istmpfilepreserved"]
			rm(tmpdir, recursive=true, force=true)
		end

    elseif  downloadtype == "Earthquake" || downloadtype == "earthquake"

		@warn("temporarily deprecated and will be implemented using SeisIO.Quake module.")

    else
		#@println("Download type is not known (chose Noise or Earthquake).")
		@error("Download type is not known (only available `Noise` for the moment).")
    end

    return 0

end
