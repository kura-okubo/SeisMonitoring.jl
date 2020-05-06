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

    if RawData_path = "default"
        RawData_path = joinpath(InputDict["fodir"], "RawData.jld2")
    else
        RawData_path = InputDict["RawData_path"]
    end

    println("***************************************")
    println("IsKurtosisRemoval  = $(InputDict["IsKurtosisRemoval"])"
    println("IsSTALTARemoval    = $(InputDict["IsSTALTARemoval"])")
    println("IsWhitening        = $(InputDict["IsWhitening"])")
    println("***************************************\n")

    t = jldopen(RawData_path, "r")

    if !haskey(t, "Waveforms")
        error("$(RawData_path) does not have Waveforms group. Please check format.")
    end

    # parallelize with keys in Rawdata.jld2 i.e. stations
    println("-------START Removing EQ--------")

    t_removeeq = @elapsed bt_time = pmap(
        x -> map_removeEQ(x, t, InputDict),
        keys(t["Waveforms"]),
    )

    JLD2.close(t)

    println("-------START Converting--------")

    t_convert = @elapsed convert_tmpfile(InputDict)

    mean_kurtosis_cputime = mean(x->x[1]).(bt_time)
    mean_stalta_cputime   = mean(x->x[2]).(bt_time)

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("time to mean kurtosis cputime  =$(mean_kurtosis_cputime)[s]")
    println("time to mean stalta cputime    =$(mean_stalta_cputime)[s]")
    println("time to remove EQ total    =$(t_removeeq)[s]")
    println("time to convert total      =$(t_convert)[s]")

end
