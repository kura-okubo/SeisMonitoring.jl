include("map_seisxcorrelation.jl")
include("seisxcorrelation_utils.jl")
include("append_wtcorr.jl")
include("assemble_seisdata.jl")

"""
    seisxcorrelation(InputDict::Dict)

Compute cross-correlation save data in jld2 file with CorrData format.

# Arguments
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters
See example to make InputDict.
"""
function seisxcorrelation(InputDict_origin::OrderedDict)

    InputDict = parse_inputdict(InputDict_origin)

    project_outputdir = abspath(InputDict["project_outputdir"])
    InputDict["fodir"] = joinpath(project_outputdir, "cc")
    tmpdir = joinpath(project_outputdir, "cc", "cc_tmp")
    InputDict["tmpdir_cc"] = tmpdir

    if ispath(tmpdir)
        rm(tmpdir, recursive = true)
    end
    mkdir(tmpdir)

    if lowercase(InputDict["cc_RawData_path"]) == "default"
        InputDict["cc_absolute_RawData_path"] = joinpath(project_outputdir, "seismicdata", "EQRemovedData.jld2")
    else
        InputDict["cc_absolute_RawData_path"] = abspath(InputDict["cc_RawData_path"])
    end

    println("***************************************")
    println("Cross-correlation starttime= $(InputDict["starttime"])")
    println("Cross-correlation endtime = $(InputDict["endtime"])")
    println("Cross-correlation method = $(InputDict["cc_method"])")
    println("Station pairs option = $(InputDict["pairs_option"])")
    println("***************************************\n")

    # get cc time windows
    InputDict["starts"], InputDict["ends"] = get_cc_time_windows(InputDict["cc_time_unit"], InputDict["sampling_frequency"],
                                        InputDict["starttime"], InputDict["endtime"])
    # scan station info
    StationDict = scan_stations(InputDict["cc_absolute_RawData_path"])
    # get station pairs
    StationPairDict = get_stationpairs(StationDict, InputDict["cc_method"], InputDict["pairs_option"])

    println("-------START Cross-correlation--------")

    # Parallelize with stationpairs
    t_removeeq = @elapsed pmap(x->map_seisxcorrelation(x, StationPairDict, InputDict), keys(StationPairDict))

    println("seisxcorrelation has been successfully done.")
    rm(tmpdir, recursive=true, force=true)

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("time for cross-correlation  =$(t_removeeq)[s]\n")


end
