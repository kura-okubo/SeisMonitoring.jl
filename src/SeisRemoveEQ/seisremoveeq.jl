include("map_removeEQ.jl")

"""
	seisremoveeq(InputDict::OrderedDict)

remove transient signals using STA/LTA and Kurtosis
"""
function seisremoveeq(InputDict_origin::OrderedDict)

	# #===================================================================#
	# #NOTE: Exit if nprocs is larger than 64 to avoid following issues:
	# # 1. massive HTTP request causing degradation of data server.
	# # 2. massive File I/O
	# Updated: this is resolved by saving seischannel files in SeisIO native format.
	# #===================================================================#
	#
	# if nprocs() > 64
	# 	error("nprocs $(nprocs()) should be less than 64 to avoid massive FIle I/O.")
	# end

    InputDict = parse_inputdict(InputDict_origin)

    project_outputdir = abspath(InputDict["project_outputdir"])
    InputDict["fodir"] = joinpath(project_outputdir, "seismicdata")
	# tmpdir = joinpath(project_outputdir, "seismicdata", "seisremoveeq_tmp")
	tmpdir = joinpath(project_outputdir, "seismicdata", "seisremoveeq")
    InputDict["tmpdir"] = tmpdir

    # if ispath(tmpdir)
    #     rm(tmpdir, recursive = true)
    # end
    !ispath(tmpdir) && mkdir(tmpdir)

    if lowercase(InputDict["RawData_path"]) == "default"
        InputDict["RawData_path"] = joinpath(InputDict["fodir"], "rawseismicdata")
    end

    println("***************************************")
    println("IsKurtosisRemoval  = $(InputDict["IsKurtosisRemoval"])")
    println("IsSTALTARemoval    = $(InputDict["IsSTALTARemoval"])")
    println("IsWhitening        = $(InputDict["IsWhitening"])")
    println("***************************************\n")

	# DEBUG: fi cannot be passed to pmap function due to parallel read issue. please reload in the pmap function
    # ispath(InputDict["RawData_path"]) ? (fi = jldopen(InputDict["RawData_path"], "r")) : error("$(InputDict["RawData_path"]) is not found.")
    # !haskey(fi, "Waveforms") && error("$(InputDict["RawData_path"]) does not have Waveforms group. Please check the waveform data format in JLD2.")
	# stations = keys(fi["Waveforms"])
 	# JLD2.close(fi)

    # parallelize with each seismic data files
	# rawdata_path_all = SeisIO.ls(InputDict["RawData_path"])
	# NOTE: updated for hierarchical directory tree 2020.10.04
	rawdata_path_all = []
	for (root, dirs, files) in ScanDir.walkdir(InputDict["RawData_path"])
       for file in files
		   fi = joinpath(root, file)
		   (split(fi, ".")[end] == "seisio") && push!(rawdata_path_all, fi)# filter if it is .seisio
       end
    end
	# println(rawdata_path_all)

    println("-------START Removing EQ--------")

    t_removeeq = @elapsed bt_time = pmap(x -> map_removeEQ(x, InputDict), rawdata_path_all)

    println("-------START Converting--------")

	# t_convert = @elapsed convert_tmpfile_seisremoveeq(InputDict)
	# t_convert = @elapsed convert_tmpfile(InputDict, "seisremoveeq")
	t_convert = 0.0
    mean_kurtosis_cputime = mean((x->x[1]).(bt_time))
    mean_stalta_cputime   = mean((x->x[2]).(bt_time))

	# rm(tmpdir, recursive=true, force=true)

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("time for mean kurtosis cputime  =$(mean_kurtosis_cputime)[s]")
    println("time for mean stalta cputime    =$(mean_stalta_cputime)[s]")
    println("time for remove EQ total    =$(t_removeeq)[s]")
    println("time for convert total      =$(t_convert)[s]")

end
