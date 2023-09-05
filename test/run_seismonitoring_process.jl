using Test
using SeisIO, SeisNoise, SeisMonitoring, JLD2, DataStructures, Dates, CSV, DataFrames, Distributed

# This test checks if the ambient seismic noise processing is adequately done in a give environment.

project_name = "run_seismonitoring_test"
project_inputdir="./data"
project_outputdir="./data"
master_param="./data/mainparam_master.jl"
fo_mainparam = project_inputdir*"/$(project_name)_INPUT/mainparam.jl"

ifGenerateTrueFiles = false # set true to generate true files and false when you test

@testset "seismonitoring_processings" begin
@testset "init project" begin

    # using the work flow in tutorial notebook in https://github.com/kura-okubo/SeisMonitoring_Example

    isdir("./data/$(project_name)_INPUT") && rm("./data/$(project_name)_INPUT", recursive=true)
    isdir("./data/$(project_name)_OUTPUT") && rm("./data/$(project_name)_OUTPUT", recursive=true)

    # make default input and output directory
    init_project(project_name=project_name,
                 project_inputdir=project_inputdir,
                 project_outputdir=project_outputdir,
                 force=true, gui=false
                )

    @test isdir("./data/$(project_name)_INPUT")
    @test isdir("./data/$(project_name)_OUTPUT")

    # make request station list
    fipath = "./data/BP_gmap-stations_example.txt"

    locchan = Dict(
                "BP" => [("*", "SP1")]
                )

    station_fodir = project_outputdir*"/$(project_name)_OUTPUT"
    station_path = make_requeststation_fromIRISgmap(fipath, locchan=locchan, fodir=station_fodir, foname="$(project_name).jld2")

    # copy master parameter file into project directory
    ~isfile(fo_mainparam) && cp(master_param, fo_mainparam)

    # replace parameters for casestudy
    set_parameter(fo_mainparam, "project_name", project_name)
    set_parameter(fo_mainparam, "project_inputdir", project_inputdir*"/$(project_name)_INPUT")
    set_parameter(fo_mainparam, "project_outputdir", project_outputdir*"/$(project_name)_OUTPUT")
    set_parameter(fo_mainparam, "requeststation_file",  station_fodir*"/$(project_name).jld2")

    PARAM = include(fo_mainparam)
    @test PARAM["requeststation_file"][1] == "./data/run_seismonitoring_test_OUTPUT/run_seismonitoring_test.jld2"

    NP = 1 # >= 1, number of additional processors; for docker users, the maximum number is managed in the docker configuration.
    addprocs(NP)
    # You need to redefine the packages in all the processors.
    @everywhere using SeisMonitoring
    @test nprocs() == 2

end

@testset "seisdownload" begin

    # For the test, we minimize the download size of data: 3 minute, 20Hz of a single channel.
    set_parameter(fo_mainparam, "download_time_unit", "60")
    set_parameter(fo_mainparam, "download_margin", "1")
    set_parameter(fo_mainparam, "starttime", "2016-01-15T00:00:00")
    set_parameter(fo_mainparam, "endtime", "2016-01-15T00:03:00")

    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=true,
            run_seisremoveeq=false,
            run_seisxcorrelation=false,
            run_seisstack=false)

    # Check if the downloaded data is valid
    S = rseis("./data/run_seismonitoring_test_OUTPUT/seismicdata/rawseismicdata/BP.EADB.40.SP1/2016/BP.EADB.40.SP1__2016-01-15T00:00:00__2016-01-15T00:01:00__sp1.seisio")
    ifGenerateTrueFiles && CSV.write("./data/seisdl_Sdata_true.csv", Tables.table(S[1].x))
    Sx_true = CSV.read("./data/seisdl_Sdata_true.csv", DataFrame).Column1
    @test typeof(S[1]) == SeisChannel
    @test Sx_true ≈ S[1].x atol=1e-6

end

@testset "seisremoveeq" begin
    set_parameter(fo_mainparam, "shorttime_window", "5")
    set_parameter(fo_mainparam, "longtime_window", "60")
    set_parameter(fo_mainparam, "timewindow_overlap", "2.5")
    set_parameter(fo_mainparam, "kurtosis_threshold", "1.0")
    set_parameter(fo_mainparam, "stalta_threshold", "1.5")
    set_parameter(fo_mainparam, "fixed_tukey_margin", "2.0")
    set_parameter(fo_mainparam, "IsWhitening", "true")

    SeisMonitoring.run_job(fo_mainparam,
        run_seisdownload=false,
        run_seisremoveeq=true,
        run_seisxcorrelation=false,
        run_seisstack=false)

    Sremeq = rseis("./data/run_seismonitoring_test_OUTPUT/seismicdata/seisremoveeq/BP.EADB.40.SP1/2016/BP.EADB.40.SP1__2016-01-15T00:00:00__2016-01-15T00:01:00__sp1.seisio")
    ifGenerateTrueFiles && CSV.write("./data/seisremoveeq_Sdata_true.csv", Tables.table(Sremeq[1].x))
    Sremeq_true = CSV.read("./data/seisremoveeq_Sdata_true.csv", DataFrame).Column1
    @test typeof(Sremeq[1]) == SeisChannel
    @test Sremeq_true ≈ Sremeq[1].x atol=1e-6
end

@testset "seisxcorrelation" begin

    # 1. filter with butterworth
    set_parameter(fo_mainparam, "freqency_band", "0.9, 1.2") # boundaries of the frequency bands
    set_parameter(fo_mainparam, "cc_time_unit", "60")
    set_parameter(fo_mainparam, "cc_len", "30")
    set_parameter(fo_mainparam, "cc_step", "5")
    set_parameter(fo_mainparam, "maxlag", "10")
    set_parameter(fo_mainparam, "cc_bpfilt_method", "Butterworth")
    set_parameter(fo_mainparam, "IsOnebit", "true") # We mitigate the edge effect during the whitening of SeisRemoveEQ for this test.
    set_parameter(fo_mainparam, "IsPreStack", "true") # avoid the prestack of the unit time (e.g. day) to keep all the short-time (e.g. hourly) CFs.

    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=false,
            run_seisremoveeq=false,
            run_seisxcorrelation=true,
            run_seisstack=false)

    mv("./data/$(project_name)_OUTPUT/cc/BP.EADB.40.SP1-BP.EADB.40.SP1", "./data/$(project_name)_OUTPUT/cc/C1_BP.EADB.40.SP1-BP.EADB.40.SP1", force=true)

    t1 = jldopen("./data/$(project_name)_OUTPUT/cc/C1_BP.EADB.40.SP1-BP.EADB.40.SP1/BP.EADB.40.SP1-BP.EADB.40.SP1__2016-01-15T00:00:00__2016-01-15T00:03:00.jld2", "r")
    C1 = t1["2016-01-15T00:00:00--2016-01-15T00:01:00/0.9-1.2"]
    close(t1)
    ifGenerateTrueFiles && CSV.write("./data/seisxcorr_C1_true.csv", Tables.table(C1.corr))
    C1_true = CSV.read("./data/seisxcorr_C1_true.csv", DataFrame).Column1
    @test typeof(C1) == CorrData
    @test C1_true ≈ C1.corr atol=1e-5
    # note: the

    #2. filter with Wavelet transform
    set_parameter(fo_mainparam, "cc_bpfilt_method", "Wavelet")

    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=false,
            run_seisremoveeq=false,
            run_seisxcorrelation=true,
            run_seisstack=false)


    t2 = jldopen("./data/$(project_name)_OUTPUT/cc/BP.EADB.40.SP1-BP.EADB.40.SP1/BP.EADB.40.SP1-BP.EADB.40.SP1__2016-01-15T00:00:00__2016-01-15T00:03:00.jld2", "r")
    C2 = t2["2016-01-15T00:00:00--2016-01-15T00:01:00/0.9-1.2"]
    close(t2)
    ifGenerateTrueFiles && CSV.write("./data/seisxcorr_C2_true.csv", Tables.table(C2.corr))
    C2_true = CSV.read("./data/seisxcorr_C2_true.csv", DataFrame).Column1
    @test typeof(C2) == CorrData
    @test C2_true ≈ C2.corr atol=1e-5
    @test C1.corr != C2.corr

    # lagt = -C1.maxlag:1/(C1.fs):C1.maxlag
    # plot(lagt, C1.corr) # The asymmetric ACF might be due to the ButterWorth filter even we apply it with zerophase=true.
    # plot!(lagt, C2.corr)

end

@testset "seisstack" begin

    # Channel collection
    set_parameter(fo_mainparam, "collect_stationpairs", "true") # collect the stationpair
    set_parameter(fo_mainparam, "compute_reference", "true") # compute the reference stack
    set_parameter(fo_mainparam, "compute_shorttimestack", "true") # compute the short time (e.g. monthly) stack

    set_parameter(fo_mainparam, "stack_RawData_dir", "./data/$(project_name)_OUTPUT/cc") # collect the stationpair
    set_parameter(fo_mainparam, "codaslice_debugplot", "false") # debug plot of the coda window
    set_parameter(fo_mainparam, "use_local_tmpdir", "false") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.

    set_parameter(fo_mainparam, "averagestack_factor", "2") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.
    set_parameter(fo_mainparam, "averagestack_step", "1") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.

    set_parameter(fo_mainparam, "reference_starttime", "2016-01-15T00:00:00") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.
    set_parameter(fo_mainparam, "reference_endtime", "2016-01-15T00:05:00") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.

    set_parameter(fo_mainparam, "min_ballistic_twin", "2") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.
    set_parameter(fo_mainparam, "max_coda_length", "8") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.
    set_parameter(fo_mainparam, "keep_corrtrace", "true") # used mainly when you run in the cluster; process the data after copying from scratch filesystem to the local machine to avoid frequent I/O to the scratch.

    # Test MWCS
    set_parameter(fo_mainparam, "measurement_method", "mwcs") # select the method of the measurement of dv/v
    !isdir("./data/$(project_name)_OUTPUT/stack") && mkdir("./data/$(project_name)_OUTPUT/stack")
    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=false,
            run_seisremoveeq=false,
            run_seisxcorrelation=false,
            run_seisstack=true)

    mv("./data/$(project_name)_OUTPUT/stack", "./data/$(project_name)_OUTPUT/stack_MWCS", force=true)

    t3 = jldopen("./data/$(project_name)_OUTPUT/stack_MWCS/reference/reference_BP.EADB-BP.EADB-11.jld2", "r")
    Cref = t3["2016-01-15T00:00:00--2016-01-15T00:05:00/0.9-1.2"]
    ifGenerateTrueFiles && CSV.write("./data/seisstack_Cref_true.csv", Tables.table(Cref.corr))
    Cref_true = CSV.read("./data/seisstack_Cref_true.csv", DataFrame).Column1
    @test Cref_true ≈ Cref.corr atol=1e-4

    t4 = jldopen("./data/$(project_name)_OUTPUT/stack_MWCS/shorttime/shorttime_BP.EADB-BP.EADB-11.jld2", "r")
    Css = t4["2016-01-15T00:00:00--2016-01-15T00:02:00/0.9-1.2"]
    ifGenerateTrueFiles && CSV.write("./data/seisstack_Css_MWCS_true.csv", Tables.table(Css.corr))
    Css_true = CSV.read("./data/seisstack_Css_MWCS_true.csv", DataFrame).Column1
    @test Css_true ≈ Css.corr atol=1e-4

    # lagt = -Css.maxlag:1/(Css.fs):Css.maxlag
    # plot(lagt, Cref.corr) # The asymmetric ACF might be due to the ButterWorth filter even we apply it with zerophase=true.
    # plot!(lagt, Css.corr)

    @test Css.misc["dvv_mwcs"] ≈ 0.01028110015 atol=1e-5
    @test Css.misc["dvv0_mwcs"] ≈ 0.01019001658 atol=1e-5
    @test Css.misc["dvv_err_mwcs"] ≈ 5.7833216142e-6 atol=1e-7

    # Test Stretching
    set_parameter(fo_mainparam, "measurement_method", "stretching") # select the method of the measurement of dv/v
    !isdir("./data/$(project_name)_OUTPUT/stack") && mkdir("./data/$(project_name)_OUTPUT/stack")
    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=false,
            run_seisremoveeq=false,
            run_seisxcorrelation=false,
            run_seisstack=true)

    mv("./data/$(project_name)_OUTPUT/stack", "./data/$(project_name)_OUTPUT/stack_Stretching", force=true)

    t3 = jldopen("./data/$(project_name)_OUTPUT/stack_Stretching/reference/reference_BP.EADB-BP.EADB-11.jld2", "r")
    Cref = t3["2016-01-15T00:00:00--2016-01-15T00:05:00/0.9-1.2"]
    ifGenerateTrueFiles && CSV.write("./data/seisstack_Cref_true.csv", Tables.table(Cref.corr))
    Cref_true = CSV.read("./data/seisstack_Cref_true.csv", DataFrame).Column1
    @test Cref_true ≈ Cref.corr atol=1e-4

    t4 = jldopen("./data/$(project_name)_OUTPUT/stack_Stretching/shorttime/shorttime_BP.EADB-BP.EADB-11.jld2", "r")
    Css = t4["2016-01-15T00:00:00--2016-01-15T00:02:00/0.9-1.2"]
    ifGenerateTrueFiles && CSV.write("./data/seisstack_Css_Stretching_true.csv", Tables.table(Css.corr))
    Css_true = CSV.read("./data/seisstack_Css_Stretching_true.csv", DataFrame).Column1
    @test Css_true ≈ Css.corr atol=1e-4

    # lagt = -Css.maxlag:1/(Css.fs):Css.maxlag
    # plot(lagt, Cref.corr) # The asymmetric ACF might be due to the ButterWorth filter even we apply it with zerophase=true.
    # plot!(lagt, Css.corr)

    @test Css.misc["dvv_ts"] ≈ -1.57503006 atol=1e-2
    @test Css.misc["cc_ts"] ≈ 0.82652899 atol=1e-5
    @test Css.misc["err_ts"] ≈ 1.13953928 atol=1e-5

end

@testset "seisstats" begin
    # Output the CSV file of the dv/v measurement
    fodir = "./data"
    starttime=DateTime("2016-1-15")
    endtime=DateTime("2016-1-16")

    # Use smstats_read_stretching or smstats_read_mwcs corresponding to the measurement method of dv/v
    SeisMonitoring.smstats_read_mwcs("./data/$(project_name)_OUTPUT/stack_MWCS/shorttime", fodir, starttime, endtime, foname="run_seismonitoring_test_MWCS.csv")
    SeisMonitoring.smstats_read_stretching("./data/$(project_name)_OUTPUT/stack_Stretching/shorttime", fodir, starttime, endtime, foname="run_seismonitoring_test_Stretching.csv")

    ifGenerateTrueFiles && cp("./data/run_seismonitoring_test_MWCS.csv", "./data/run_seismonitoring_test_MWCS_true.csv", force=true)
    ifGenerateTrueFiles && cp("./data/run_seismonitoring_test_Stretching.csv", "./data/run_seismonitoring_test_Stretching_true.csv", force=true)
    df_dvv_MWCS = CSV.read("./data/run_seismonitoring_test_MWCS.csv", DataFrame)
    df_dvv_MWCS_true = CSV.read("./data/run_seismonitoring_test_MWCS_true.csv", DataFrame)
    @test df_dvv_MWCS_true[!, "dvv_mwcs"] ≈ df_dvv_MWCS_true[!, "dvv_mwcs"] atol=1e-7
    @test df_dvv_MWCS_true[!, "dvv0_mwcs"] ≈ df_dvv_MWCS_true[!, "dvv0_mwcs"] atol=1e-7
    @test df_dvv_MWCS_true[!, "dvv_err_mwcs"] ≈ df_dvv_MWCS_true[!, "dvv_err_mwcs"] atol=1e-9

    df_dvv_Stretching = CSV.read("./data/run_seismonitoring_test_Stretching.csv", DataFrame)
    df_dvv_Stretching_true = CSV.read("./data/run_seismonitoring_test_Stretching_true.csv", DataFrame)
    @test df_dvv_Stretching_true[!, "dvv_ts"] ≈ df_dvv_Stretching[!, "dvv_ts"] atol=1e-6
    @test df_dvv_Stretching_true[!, "cc_ts"] ≈ df_dvv_Stretching[!, "cc_ts"] atol=1e-6
    @test df_dvv_Stretching_true[!, "err_ts"] ≈ df_dvv_Stretching[!, "err_ts"] atol=1e-5

end
end
