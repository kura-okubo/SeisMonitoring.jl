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
    # StationDict, all_stations = scan_stations(InputDict["cc_absolute_RawData_path"])
    # println(StationDict)
    # println(all_stations)
    # get station pairs
    # netstachan1_list, netstachan2_list, StationPairs = get_stationpairs(StationDict, InputDict["cc_normalization"], InputDict["pairs_option"], InputDict["chanpair_type"])

    println("-------START Cross-correlation--------")

    # compute time chunk to be parallelized

    t_assemble_all = []
    t_fft_all = []
    t_corr_all = []

    rawdata_path_all = SeisIO.ls(InputDict["cc_absolute_RawData_path"])

    # output timechunk cpu time
    fi_chunkcpu = Base.open(joinpath(project_outputdir, "timechunk_cputime.txt"), "w")

    for timechunkid in Iterators.partition(1:length(InputDict["starts"]), InputDict["timechunk_increment"])
    # for timechunkid in Iterators.partition(1:length(InputDict["starts"]), 10)

        ct_1 = now()

        InputDict["starts_chunk"] = InputDict["starts"][timechunkid]
        InputDict["ends_chunk"] = InputDict["ends"][timechunkid]
        println("$(now()): time chunk $(u2d(InputDict["starts_chunk"][1]))-$(u2d(InputDict["ends_chunk"][end]))")
        # 1. compute FFT and store data into FFTDict

        # make chunk_fi_stationdict to scan stations
        InputDict["chunk_fi_stationdict"] = get_chunk_fi_stationdict(rawdata_path_all,
                                                u2d(InputDict["starts_chunk"][1]), u2d(InputDict["ends_chunk"][end]))

        all_stations_chunk = sort(collect(keys(InputDict["chunk_fi_stationdict"]))) # list of all station channels within this timechunk.
        # scan
        netstachan1_list, netstachan2_list, StationPairs_chunk = get_stationpairs_chunk(all_stations_chunk,
                                                             InputDict["cc_normalization"], InputDict["pairs_option"], InputDict["chanpair_type"])

        @show all_stations_chunk
        @show length(StationPairs_chunk)
        #NOTE:
        # Specify WorkerPool to avoid the following issue on redundant precompile through chunk loop:
        # e.g. Given 100 cores as workers, and parallelize 10 stations for fft:
        #
        # - first chunk: worker 1,2 ... 10 precompiled
        # - second chunk: worker 1,2, ..., 9, 11, taking precompile time for worker 11
        # - third chunk: worker 1,2, ..., 9, 12, taking precompile time for worker 12
        # - fourth chunk: worker 1,2, ..., 9, 13, taking precompile time for worker 13

        N_workpool_fft = min(length(all_stations_chunk)+1, nworkers()+1)
        map_compute_fft_workerpool = WorkerPool(collect(2:N_workpool_fft))

        N_workpool_cc = min(length(StationPairs_chunk)+1, nworkers()+1)
        map_compute_cc_workerpool = WorkerPool(collect(2:N_workpool_cc))


        let FFTs, FFT_Dict

            ta_1 = @elapsed A = pmap(x -> map_compute_fft(x, InputDict), map_compute_fft_workerpool, all_stations_chunk) # store FFTs in memory and deallocate after map_compute_correlation().
            stations    = (x->x[1]).(A)
            FFTs        = (x->x[2]).(A)
            push!(t_assemble_all, mean((x->x[3]).(A)))
            push!(t_fft_all, mean((x->x[4]).(A)))

            FFT_Dict = Dict{String,Dict{String, FFTData}}()

            for (i, station) in enumerate(stations)
                FFT_Dict[station] = FFTs[i]
            end

            # memory_use=sizeof(FFTs)/1e9 #[GB]
            memory_use=Base.summarysize(FFTs)/1e9
            println("debug: memory_use: $(memory_use) GB.")
            memory_use > InputDict["MAX_MEM_USE"] && @error("Memory use during FFT exceeds MAX_MEM_USE ($(memory_use)GB is used). Please decrease timechunk_increment.")

            #NOTE: using map() function to move each pair of FFTData from host to workers.
            # ta_2 = @elapsed B = pmap((x, y) -> map_compute_cc(x, y, InputDict),
            #                                 map_compute_cc_workerpool,
            #                                 map((k, l) -> (FFT_Dict[k], FFT_Dict[l]), netstachan1_list, netstachan2_list),
            #                                     StationPairs_chunk)

            tm1 = @elapsed FFT1_dict = map(k -> FFT_Dict[k], netstachan1_list)
            tm2 = @elapsed FFT2_dict = map(l -> FFT_Dict[l], netstachan2_list)
            println("$(now()): tmap1, tmap2 = $(tm1), $(tm2)[s]")

            # ta_2 = @elapsed B = pmap((fft1, fft2, pair) -> map_compute_cc(fft1, fft2, pair, InputDict),
            #                                 map_compute_cc_workerpool,
            #                                 FFT1_dict, FFT2_dict, StationPairs_chunk)
            ta_2 = 0 #DEBUG
            ta_3 = @elapsed C = pmap(x -> pmaptest_1(x), map_compute_cc_workerpool, StationPairs_chunk)

            # push!(t_corr_all, mean((x->x[1]).(B)))
            push!(t_corr_all, 0) #DEBUG
            println("time for map_fft, map_cc, testpmap = $(ta_1), $(ta_2), $(ta_3) [s]")
        end

        ct_2 = now()
        ct_elapse = ct_2-ct_1

        Base.write(fi_chunkcpu, "$(ct_elapse.value)\n")

    end

    Base.close(fi_chunkcpu)
    # Parallelize with stationpairs
    # t_removeeq = @elapsed bt_time = pmap(x->map_seisxcorrelation(x, InputDict),
    #                                             StationPairs)
    #
    println("seisxcorrelation has been successfully done.")

    # # DEBUG: to avoid error in return in cluseter, comment out for the moment
    # mean_assemble_cputime   = mean((x->x[1]).(bt_time))
    # mean_fft_cputime        = mean((x->x[2]).(bt_time))
    # mean_xcorr_cputime      = mean((x->x[3]).(bt_time))
    !isempty(t_assemble_all) ? (median_assemble_all = median(t_assemble_all)) : (median_assemble_all=0)
    !isempty(t_fft_all)      ? (median_fft_all = median(t_fft_all)) : (median_fft_all=0)
    !isempty(t_corr_all)     ? (median_corr_all = median(t_corr_all)) : (median_corr_all=0)


    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("median time for assemble cputime  =$(median_assemble_all)[s]")
    println("median time for fft cputime       =$(median_fft_all)[s]")
    println("median time for cc  cputime       =$(median_corr_all)[s]")
    # println("time for mean assemble cputime     =$(mean_assemble_cputime)[s]")
    # println("time for mean fft cputime          =$(mean_fft_cputime)[s]")
    # println("time for mean cross-correlation cputime = $(mean_xcorr_cputime)[s]")
    println("done.")
end
