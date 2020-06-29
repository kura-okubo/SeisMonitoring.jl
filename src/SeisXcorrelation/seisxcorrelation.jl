include("map_seisxcorrelation.jl")
include("seisxcorrelation_utils.jl")
include("compute_frequency_decomposition.jl")
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

    if lowercase(InputDict["cc_RawData_path"]) == "default"
        InputDict["cc_absolute_RawData_path"] = joinpath(project_outputdir, "seismicdata", "EQRemovedData.jld2")
    else
        InputDict["cc_absolute_RawData_path"] = abspath(InputDict["cc_RawData_path"])
    end

    println("***************************************")
    println("Cross-correlation starttime     = $(InputDict["starttime"])")
    println("Cross-correlation endtime       = $(InputDict["endtime"])")
    println("Cross-correlation normalization = $(InputDict["cc_normalization"])")
    println("Cross-correlation type          = $(InputDict["corr_type"])")
    println("Station pairs option            = $(InputDict["pairs_option"])")
    println("Cross-correlation chanpair type = $(InputDict["chanpair_type"])")
    println("***************************************\n")

    # get cc time windows
    InputDict["starts"], InputDict["ends"] = get_cc_time_windows(InputDict["cc_time_unit"], InputDict["sampling_frequency"],
                                        InputDict["starttime"], InputDict["endtime"])
    # scan station info
    StationDict = scan_stations(InputDict["cc_absolute_RawData_path"])
    # get station pairs
    StationPairDict = get_stationpairs(StationDict, InputDict["cc_normalization"], InputDict["pairs_option"], InputDict["chanpair_type"])

    println("-------START Cross-correlation--------")

    # Parallelize with stationpairs
    t_removeeq = @elapsed bt_time = pmap(x->map_seisxcorrelation(x, StationPairDict, InputDict),
                                                    collect(keys(StationPairDict)))

    println("seisxcorrelation has been successfully done.")

    # DEBUG: to avoid error in return in cluseter, comment out for the moment
    mean_assemble_cputime   = mean((x->x[1]).(bt_time))
    mean_fft_cputime        = mean((x->x[2]).(bt_time))
    mean_xcorr_cputime      = mean((x->x[3]).(bt_time))

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("Total time for cross-correlation   =$(t_removeeq)[s]")
    println("time for mean assemble cputime     =$(mean_assemble_cputime)[s]")
    println("time for mean fft cputime          =$(mean_fft_cputime)[s]")
    println("time for mean cross-correlation cputime = $(mean_xcorr_cputime)[s]")

end
