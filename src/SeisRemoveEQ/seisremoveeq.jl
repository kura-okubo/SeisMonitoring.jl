include("map_removeEQ.jl")

"""
	seisremoveeq(InputDict::OrderedDict)

remove transient signals using STA/LTA and Kurtosis
"""
function seisremoveeq(InputDict_origin::OrderedDict)

    InputDict = parse_inputdict(InputDict_origin)

    project_outputdir = abspath(InputDict["project_outputdir"])
    InputDict["fodir"] = joinpath(project_outputdir, "seismicdata")
    tmpdir = joinpath(project_outputdir, "seismicdata", "seisremoveeq_tmp")
    InputDict["tmpdir_rem"] = tmpdir

    if ispath(tmpdir)
        rm(tmpdir, recursive = true)
    end
    mkdir(tmpdir)

    if lowercase(InputDict["RawData_path"]) == "default"
        InputDict["RawData_path"] = joinpath(InputDict["fodir"], "RawData.jld2")
    end

    println("***************************************")
    println("IsKurtosisRemoval  = $(InputDict["IsKurtosisRemoval"])")
    println("IsSTALTARemoval    = $(InputDict["IsSTALTARemoval"])")
    println("IsWhitening        = $(InputDict["IsWhitening"])")
    println("***************************************\n")

	#DEBUG: fi cannot be passed to pmap function due to parallel read issue. please reload in the pmap function
    ispath(InputDict["RawData_path"]) ? (fi = jldopen(InputDict["RawData_path"], "r")) : error("$(InputDict["RawData_path"]) is not found.")
    !haskey(fi, "Waveforms") && error("$(RawData_path) does not have Waveforms group. Please check the waveform data format in JLD2.")
	stations = keys(fi["Waveforms"])
 	JLD2.close(fi)

    # parallelize with keys in Rawdata.jld2 i.e. stations
    println("-------START Removing EQ--------")

    t_removeeq = @elapsed bt_time = pmap(x -> map_removeEQ(x, InputDict),stations)

    JLD2.close(t)

    println("-------START Converting--------")

    t_convert = @elapsed convert_tmpfile_seisremoveeq(InputDict)

    mean_kurtosis_cputime = mean((x->x[1]).(bt_time))
    mean_stalta_cputime   = mean((x->x[2]).(bt_time))

	rm(tmpdir, recursive=true, force=true)

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("time for mean kurtosis cputime  =$(mean_kurtosis_cputime)[s]")
    println("time for mean stalta cputime    =$(mean_stalta_cputime)[s]")
    println("time for remove EQ total    =$(t_removeeq)[s]")
    println("time for convert total      =$(t_convert)[s]")

end
