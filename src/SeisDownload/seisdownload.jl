include("seisdownload_utils.jl")
include("map_seisdownload.jl")

"""
    seisdownload(InputDict::Dict)

    Request seismic data and save into jld2.
# Arguments
- `InputDict`    : dictionary which contains request information
"""
function seisdownload(InputDict_origin::OrderedDict)

	#===================================================================#
	#NOTE: Exit if nprocs is larger than 64 to avoid following issues:
	# 1. massive HTTP request causing degradation of data server.
	# 2. massive File I/O
	#===================================================================#

	if nprocs() > 64
		error("nprocs $(nprocs()) should be less than 64 to avoid massive HTTP request and FIle I/O.")
	end

	# parse input dictionary
	InputDict = parse_inputdict(InputDict_origin)

	project_outputdir		= abspath(InputDict["project_outputdir"])
	fodir					= joinpath(project_outputdir, "seismicdata")
	tmpdir 					= joinpath(project_outputdir, "seismicdata", "rawseismicdata")
	InputDict["fodir"] 		= fodir
	InputDict["tmpdir"] 	= tmpdir

	stationxml_dir = joinpath(fodir, "stationxml")
	InputDict["stationxml_dir"] = stationxml_dir


	# if ispath(tmpdir); rm(tmpdir, recursive=true); end
	#NOTE: do not delete output file dir for restarting job
	!ispath(tmpdir) && mkdir(tmpdir)

	# if ispath(stationxml_dir); rm(stationxml_dir, recursive=true); end
	!ispath(stationxml_dir) && mkdir(stationxml_dir)

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

	# calculate start time list (starttimelist) with each Download_time_unit
	starttimelist = get_starttimelist(starttime, endtime, download_time_unit)
	# generate DLtimestamplist and ststationlist
	# DLtimestamplist = Utils.get_timestamplist(starttimelist)

	InputDict["starttimelist"] = starttimelist
	# InputDict["DLtimestamplist"] = DLtimestamplist

	#----Restrict number of processors------#
	#NEVER CHANGE THIS THRESHOLD OTHERWISE IT OVERLOADS THE DATA SERVER
	np = nprocs()
	if np > 64 throw(DomainError(np, "np must be smaller than 64.")) end
	#---------------------------------------#

    # # Test download to evaluate use of memory and estimate download time.
	# InputDict_test = deepcopy(InputDict) # to avoid overwriting InputDict; unknown bug while deepcopying in testdownload function
	# testdownload(InputDict_test, length(starttimelist))

	# Start downloading data
	println("-------START DOWNLOADING--------")
	t_download = @elapsed pmap(x -> map_seisdownload_NOISE(x, InputDict), 1:length(starttimelist))

	# convert intermediate file to prescibed file format (JLD2, ASDF, ...)

	println("-------START CONVERTING--------")

	# t_convert = @elapsed convert_tmpfile(InputDict)
	# t_convert = @elapsed convert_tmpfile(InputDict, "seisdownload")


	println("---Summary of computational time---")
	println(@sprintf("Total download time:%8.4f[s]", t_download))
	# println(@sprintf("Total convert time:%8.4f[s]", t_convert))

	# if !InputDict["Istmpfilepreserved"]
	# 	rm(tmpdir, recursive=true, force=true)
	# end

	if !InputDict["IsXMLfilepreserved"]
		rm(stationxml_dir, recursive=true, force=true)
	end

    return 0

end
