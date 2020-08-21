# include("map_seisxcorrelation.jl")
include("map_compute_fft.jl")
include("map_compute_cc.jl")
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
        # InputDict["cc_absolute_RawData_path"] = joinpath(project_outputdir, "seismicdata", "EQRemovedData.jld2")
        InputDict["cc_absolute_RawData_path"] = joinpath(project_outputdir, "seismicdata", "seisremoveeq")
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
    StationDict, all_stations = scan_stations(InputDict["cc_absolute_RawData_path"])
    println(StationDict)
    println(all_stations)
    # get station pairs
    netstachan1_list, netstachan2_list, StationPairs = get_stationpairs(StationDict, InputDict["cc_normalization"], InputDict["pairs_option"], InputDict["chanpair_type"])

    println("-------START Cross-correlation--------")

    # compute time chunk to be parallelized
    f_debug(x, y) = println((x, y))
    function f_debug2(x::Tuple{Dict{String,SeisNoise.FFTData},Dict{String,SeisNoise.FFTData}}, y::String)
        print("debug:$(y) ")
        FFT_1, FFT_2 = x
        (isempty(FFT_1) || isempty(FFT_2)) && return
        FFT_1_1 = FFT_1[collect(keys(FFT_1))[1]]
        FFT_2_1 = FFT_2[collect(keys(FFT_2))[1]]
        id1 = FFT_1_1.name*"-"*FFT_1_1.id
        id2 = FFT_2_1.name*"-"*FFT_2_1.id
        println("should be $(id1)__$(id2)")
    end

    t_assemble_all, t_fft_all, t_corr_all = zeros(3)

    # for timechunkid in Iterators.partition(1:length(InputDict["starts"]), InputDict["timechunk_increment"])
    for timechunkid in Iterators.partition(1:length(InputDict["starts"]), 10)

        InputDict["starts_chunk"] = InputDict["starts"][timechunkid]
        InputDict["ends_chunk"] = InputDict["ends"][timechunkid]
        println("time chunk $(u2d(InputDict["starts_chunk"][1]))-$(u2d(InputDict["ends_chunk"][end]))")
        # 1. compute FFT and store data into FFTDict
        let FFTs, FFT_Dict

            A = pmap(x -> map_compute_fft(x, InputDict), all_stations) # store FFTs in memory and deallocate after map_compute_correlation().
            stations    = (x->x[1]).(A)
            FFTs        = (x->x[2]).(A)
            t_assemble_all  += sum((x->x[3]).(A))
            t_fft_all       += sum((x->x[4]).(A))
            #
            # println(stations)
            # println(typeof(FFTs))

            FFT_Dict = Dict{String,Dict}()

            for (i, station) in enumerate(stations)
                FFT_Dict[station] = FFTs[i]
            end

            #debug
            # FFT_1 = FFT_Dict[collect(keys(FFT_Dict))[2]]
            # FFT_2 = FFT_1[collect(keys(FFT_1))[1]]
            # println(FFT_2)

            # memory_use=sizeof(FFTs)/1e9 #[GB]
            memory_use=Base.summarysize(FFTs)/1e9
            println("debug: memory_use: $(memory_use) GB.")
            memory_use > InputDict["MAX_MEM_USE"] && @error("Memory use during FFT exceeds MAX_MEM_USE ($(memory_use)GB is used). Please decrease timechunk_increment.")
            # pmap((x, y) -> f_debug2(x, y),map((x, y) -> (FFT_Dict[x], FFT_Dict[y]), netstachan1_list, netstachan2_list),
            #                                 StationPairs)
            #NOTE: using map() function to move each pair of FFTData from host to workers.
            B = pmap((x, y) -> map_compute_cc(x, y, InputDict),
                                            map((x, y) -> (FFT_Dict[x], FFT_Dict[y]), netstachan1_list, netstachan2_list),
                                            StationPairs)
            t_corr_all += sum((x->x[1]).(B))

        end

    end
    # Parallelize with stationpairs
    # t_removeeq = @elapsed bt_time = pmap(x->map_seisxcorrelation(x, InputDict),
    #                                             StationPairs)
    #
    println("seisxcorrelation has been successfully done.")

    # # DEBUG: to avoid error in return in cluseter, comment out for the moment
    # mean_assemble_cputime   = mean((x->x[1]).(bt_time))
    # mean_fft_cputime        = mean((x->x[2]).(bt_time))
    # mean_xcorr_cputime      = mean((x->x[3]).(bt_time))

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("Total time for assemble cputime     =$(t_assemble_all)[s]")
    println("Total time for fft cputime          =$(t_fft_all)[s]")
    println("Total time for cross-correlation cputime = $(t_corr_all)[s]")
    # println("time for mean assemble cputime     =$(mean_assemble_cputime)[s]")
    # println("time for mean fft cputime          =$(mean_fft_cputime)[s]")
    # println("time for mean cross-correlation cputime = $(mean_xcorr_cputime)[s]")
    println("done.")
end
