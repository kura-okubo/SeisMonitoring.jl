using Test
using SeisIO, SeisNoise, JLD2, DataStructures, Dates, CSV, DataFrames, Distributed
using SeisMonitoring

# auxiliary tests
project_name = "run_seismonitoring_test"
project_inputdir="./data"
project_outputdir="./data"
master_param="./data/mainparam_master.jl"
fo_mainparam = project_inputdir*"/$(project_name)_INPUT/mainparam.jl"


ifGenerateTrueFiles = false # set true to generate true files and false when you test

@testset "dataavailability" begin

    datatype = "rawseismicdata" # choose "rawseismicdata" or "seisremoveeq"

    fidir = "./data/run_seismonitoring_test_OUTPUT/seismicdata/"*datatype
    foname = "./data/dataavailability_"*datatype

    starttime = DateTime(2016, 1, 15)
    endtime   = DateTime(2016, 1, 16)

    smstats_dataavailability(fidir,foname, starttime, endtime)

    ifGenerateTrueFiles && cp("./data/dataavailability_rawseismicdata.csv", "./data/dataavailability_rawseismicdata_true.csv", force=true)
    df_avail = CSV.read("./data/dataavailability_rawseismicdata.csv", DataFrame).data_fraction
    df_avail_true = CSV.read("./data/dataavailability_rawseismicdata_true.csv", DataFrame).data_fraction
    @test df_avail == df_avail_true
end

@testset "seisremoveeq append traces" begin

    set_parameter(fo_mainparam, "Append_alltraces", "true")

    SeisMonitoring.run_job(fo_mainparam,
        run_seisdownload=false,
        run_seisremoveeq=true,
        run_seisxcorrelation=false,
        run_seisstack=false)

    tremall = jldopen("./data/run_seismonitoring_test_OUTPUT/seismicdata/seisremoveeq/BP.EADB.40.SP1__2016-01-15T00:00:00__2016-01-15T00:01:00__sp1.seisio.jld2", "r")
    S = tremall["S"]
    @test S.misc["data_fraction"] ≈ 0.874271 atol=1e-3
    close(tremall)

end

@testset "CC spectral normalization" begin

    # Coherence method
    set_parameter(fo_mainparam, "cc_normalization", "coherence")

    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=false,
            run_seisremoveeq=false,
            run_seisxcorrelation=true,
            run_seisstack=false)

    mv("./data/$(project_name)_OUTPUT/cc/BP.EADB.40.SP1-BP.EADB.40.SP1", "./data/$(project_name)_OUTPUT/cc/Coh_BP.EADB.40.SP1-BP.EADB.40.SP1", force=true)

    t1 = jldopen("./data/$(project_name)_OUTPUT/cc/Coh_BP.EADB.40.SP1-BP.EADB.40.SP1/BP.EADB.40.SP1-BP.EADB.40.SP1__2016-01-15T00:00:00__2016-01-15T00:03:00.jld2", "r")
    C1 = t1["2016-01-15T00:00:00--2016-01-15T00:01:00/0.9-1.2"]
    close(t1)
    ifGenerateTrueFiles && CSV.write("./data/seisxcorr_Coh_true.csv", Tables.table(C1.corr))
    C1_true = CSV.read("./data/seisxcorr_Coh_true.csv", DataFrame).Column1
    @test typeof(C1) == CorrData
    @test C1_true ≈ C1.corr atol=1e-2

    # Deconvolution method
    set_parameter(fo_mainparam, "cc_normalization", "deconvolution")

    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=false,
            run_seisremoveeq=false,
            run_seisxcorrelation=true,
            run_seisstack=false)

    mv("./data/$(project_name)_OUTPUT/cc/BP.EADB.40.SP1-BP.EADB.40.SP1", "./data/$(project_name)_OUTPUT/cc/Cdec_BP.EADB.40.SP1-BP.EADB.40.SP1", force=true)

    t1 = jldopen("./data/$(project_name)_OUTPUT/cc/Cdec_BP.EADB.40.SP1-BP.EADB.40.SP1/BP.EADB.40.SP1-BP.EADB.40.SP1__2016-01-15T00:00:00__2016-01-15T00:03:00.jld2", "r")
    C1 = t1["2016-01-15T00:00:00--2016-01-15T00:01:00/0.9-1.2"]
    close(t1)
    ifGenerateTrueFiles && CSV.write("./data/seisxcorr_Cdec_true.csv", Tables.table(C1.corr))
    C1_true = CSV.read("./data/seisxcorr_Cdec_true.csv", DataFrame).Column1
    @test typeof(C1) == CorrData
    @test C1_true ≈ C1.corr atol=1e-2

end

@testset "compute_dvvdqq" begin
    @testset "seisstack_dvvdqq" begin

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

        # Test compute_dvvdqq
        set_parameter(fo_mainparam, "measurement_method", "compute_dvvdqq") # select the method of the measurement of dv/v
        set_parameter(fo_mainparam, "computedqq_smoothing_windowlength", "2.0") # select the method of the measurement of dv/v
        set_parameter(fo_mainparam, "stretch_debugplot", "true") # plot for coda Q debug figures
        set_parameter(fo_mainparam, "codaslice_debugplot", "true") # test the coda slice debug plot


        !isdir("./data/$(project_name)_OUTPUT/stack") && mkdir("./data/$(project_name)_OUTPUT/stack")
        SeisMonitoring.run_job(fo_mainparam,
                run_seisdownload=false,
                run_seisremoveeq=false,
                run_seisxcorrelation=false,
                run_seisstack=true)

        mv("./data/$(project_name)_OUTPUT/stack", "./data/$(project_name)_OUTPUT/stack_conputedvvdqq", force=true)

        tqref = jldopen("./data/$(project_name)_OUTPUT/stack_conputedvvdqq/reference/reference_BP.EADB-BP.EADB-11.jld2", "r") # everytime compute the Cref, although the process is identical over different measurment method
        Cref = tqref["2016-01-15T00:00:00--2016-01-15T00:05:00/0.9-1.2"]
        ifGenerateTrueFiles && CSV.write("./data/seisstack_Cref_true.csv", Tables.table(Cref.corr))
        Cref_true = CSV.read("./data/seisstack_Cref_true.csv", DataFrame).Column1
        @test Cref_true ≈ Cref.corr atol=1e-4

        tq = jldopen("./data/$(project_name)_OUTPUT/stack_conputedvvdqq/shorttime/shorttime_BP.EADB-BP.EADB-11.jld2", "r")
        Css = tq["2016-01-15T00:00:00--2016-01-15T00:02:00/0.9-1.2"]
        ifGenerateTrueFiles && CSV.write("./data/seisstack_Css_dvvdqq_true.csv", Tables.table(Css.corr))
        Css_true = CSV.read("./data/seisstack_Css_dvvdqq_true.csv", DataFrame).Column1
        @test Css_true ≈ Css.corr atol=1e-4


        # lagt = -Cvq.maxlag:1/(Cvq.fs):Css.maxlag
        # plot(lagt, Cref.corr) # The asymmetric ACF might be due to the ButterWorth filter even we apply it with zerophase=true.
        # plot!(lagt, Cvq.corr)

        @info("NOTE: Here we just test the consistensy in computation. The value of dqq needs to be further evaluated to check its plausibility.")
        @test Css.misc["cc_dvv"] ≈ 0.82652899 atol=1e-5
        @test Css.misc["err_dvv"] ≈ 1.1395392867 atol=1e-5
        @test Css.misc["dqq_avg"] ≈ 16.67596460789 atol=1e-2

    end

    @testset "seisstats_dvvdqq" begin
        # Output the CSV file of the dv/v measurement
        fodir = "./data"
        starttime=DateTime("2016-1-15")
        endtime=DateTime("2016-1-16")

        # Use smstats_read_dvvdqq
        smstats_read_computedvvdqq("./data/$(project_name)_OUTPUT/stack_conputedvvdqq/shorttime", fodir, starttime, endtime, foname="run_seismonitoring_test_dvvdqq.csv");

        ifGenerateTrueFiles && cp("./data/run_seismonitoring_test_dvvdqq.csv", "./data/run_seismonitoring_test_dvvdqq_true.csv", force=true)

        df_dvv_dvvdqq = CSV.read("./data/run_seismonitoring_test_dvvdqq.csv", DataFrame)
        df_dvv_dvvdqq_true = CSV.read("./data/run_seismonitoring_test_dvvdqq_true.csv", DataFrame)
        @test df_dvv_dvvdqq_true[!, "dqq_avg"] ≈ df_dvv_dvvdqq_true[!, "dqq_avg"] atol=1e-7

    end
end

@testset "stacking methods" begin

    fi = jldopen("./data/$(project_name)_OUTPUT/cc_channel_collection/BP.EADB-BP.EADB-11.jld2", "r")
    starttime=DateTime("2016-1-15")
    endtime=DateTime("2016-1-16")

    C_origin, _ = SeisMonitoring.assemble_corrdata(fi, starttime, endtime, "0.9-1.2")

    for stackmethod in ["linear", "selective", "robust", "pws", "robustpws"]

        C = deepcopy(C_origin)

        if stackmethod == "selective"
            SeisMonitoring.selectivestack(C, "shorttime")
        end

        InputDict = OrderedDict("stack_method" => stackmethod, "dist_threshold" => 1.0, "distance_type"  => "CorrDist")
        SeisMonitoring.sm_stack!(C, "shorttime", InputDict)
        @test size(C.corr, 2) == 1
        @test C.misc["stack_method"] == stackmethod

        ifGenerateTrueFiles && CSV.write("./data/seisstack_methods_$(stackmethod)_true.csv", Tables.table(C.corr))
        C_true = CSV.read("./data/seisstack_methods_$(stackmethod)_true.csv", DataFrame).Column1
        @test C_true ≈ C.corr

    end

end

@testset "seisdownload IRIS test" begin


    set_parameter(fo_mainparam, "IsLocationBox", "true")
    set_parameter(fo_mainparam, "reg",  "47.6391, 47.6462, -122.3354, -122.3229")
    set_parameter(fo_mainparam, "download_time_unit", "60")
    set_parameter(fo_mainparam, "download_margin", "1")
    set_parameter(fo_mainparam, "starttime", "2019-03-01T00:00:00")
    set_parameter(fo_mainparam, "endtime", "2019-03-01T00:01:00")
    set_parameter(fo_mainparam, "sampling_frequency", "20")
    set_parameter(fo_mainparam, "requeststation_file",  "./data/run_seismonitoring_test_OUTPUT/seismonitoring_IRIS_test.jld2")

    include("./data/$(project_name)_INPUT/mainparam.jl")

    # make station
    dftemp = DataFrame(network="UW", station="BST11",location="*",channel="HHZ")

    RequestStations = Dict(
        "IRISDMC" => dftemp
    )

    jldopen("./data/run_seismonitoring_test_OUTPUT/seismonitoring_IRIS_test.jld2", "w") do f
        for key in keys(RequestStations)
            f[key] = RequestStations[key]
        end
    end

    SeisMonitoring.run_job(fo_mainparam,
            run_seisdownload=true,
            run_seisremoveeq=false,
            run_seisxcorrelation=false,
            run_seisstack=false)

    @test isfile("./data/run_seismonitoring_test_OUTPUT/seismicdata/rawseismicdata/UW.BST11..HHZ/2019/UW.BST11..HHZ__2019-03-01T00:00:00__2019-03-01T00:00:59__hhz.seisio")

end
