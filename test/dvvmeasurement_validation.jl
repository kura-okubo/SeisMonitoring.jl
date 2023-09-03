using SeisIO, SeisNoise, JLD2, Dates, Plots, Interpolations
using SeisMonitoring: seisdvv_mwcs, seisdvv_stretching

@testset "dv/v measurement validation" begin

    # Using the work flow in
    # https://github.com/kura-okubo/SeisMonitoring_Paper/tree/develop/Others/validation_of_SeisMeasurement
    # to validate the measurement of dv/v.
    println(pwd())
    # read reference stacked trace processed in Cascadia

    t = jldopen("./data/reference_BP.EADB-BP.LCCB-11.jld2", "r")
    Cref = t["2010-01-01T00:00:00--2022-06-01T00:00:00/0.9-1.2"]
    @test typeof(Cref) == SeisNoise.CorrData

    # Workflow
    # 1. Stretch the trace from -0.1% to 0.1%
    # 2. Store it in corrdata
    # 3. Measure dv/v using seismeasurement, and see the quality of the evaluation.

    dvmin = -0.001
    dvmax = 0.001
    ntrial = 11

    tvec = -Cref.maxlag:(1/Cref.fs):Cref.maxlag
    ϵ = range(dvmin, stop=dvmax, length=ntrial)
    L = 1. .- ϵ # from velocity decrease to increase
    tau = tvec * L'

    Cstretched = deepcopy(Cref)
    Cstretched.corr = zeros(Float32, size(Cref.corr)[1], ntrial)

    # set of stretched/compressed current waveforms
    waveform_ref = Cref.corr
    for ii = 1:ntrial
    	s = LinearInterpolation(tau[:,ii],vec(waveform_ref),extrapolation_bc=Flat())(tvec)
    	Cstretched.corr[:, ii] = s
    end

    ref = vec(Cref.corr)

    # Measure dv/v using seismeasurement

    dvv_stretch = zeros(ntrial)
    dvv_mwcs = zeros(ntrial)
    dvv_mwcs0 = zeros(ntrial)

    InputDict = Dict()
    InputDict["mwcs_window_length"] = 6.0
    InputDict["mwcs_window_step"] = 3.0
    InputDict["mwcs_max_dt"] = 1.0
    InputDict["mwcs_smoothing_half_win"] = 5
    InputDict["coda_init_factor"] = 2
    InputDict["max_coda_length"] = 40
    InputDict["min_ballistic_twin"],  =  5.0
    InputDict["background_vel"] = 1000.0


    for ii = 1:ntrial
        # stretchning
        MeasurementDict_stretch = seisdvv_stretching(ref, Cstretched.corr[:, ii], Cref.misc["timelag"], Cref.misc["coda_window"],
        Cref.freqmin,Cref.freqmax,
        dvmin=-0.02,　dvmax=0.02,
        ntrial_v=101)
        dvv_stretch[ii] = 1e-2*MeasurementDict_stretch["dvv_ts"]

        # MWCS
        MeasurementDict_mwcs = seisdvv_mwcs(ref, Cstretched.corr[:, ii],Cref.freqmin,Cref.freqmax,Cref.fs,-Cref.maxlag,
        InputDict["mwcs_window_length"], InputDict["mwcs_window_step"], InputDict["mwcs_max_dt"],
        InputDict["mwcs_smoothing_half_win"], InputDict["coda_init_factor"], InputDict["max_coda_length"],
        InputDict["min_ballistic_twin"], Cref.dist*1e3, InputDict["background_vel"])

        if isempty(MeasurementDict_mwcs)
            dvv_mwcs[ii]  = 0.0
            dvv_mwcs0[ii] = 0.0
            else
            dvv_mwcs[ii]  = -MeasurementDict_mwcs["dvv_mwcs"]
            dvv_mwcs0[ii] = -MeasurementDict_mwcs["dvv0_mwcs"]
            end
    end

    # The precomputed true values
    dvv_stretch_true = [-0.0010308617234468939, -0.000823246492985972, -0.0006332665330661323, -0.00042324649298597193, -0.00019398797595190383, -2.3246492985971944e-5, 0.00016673346693386773, 0.00037675350701402805, 0.0005667334669338677, 0.000776753507014028, 0.0010060120240480961]
    dvv_mwcs_true = [-0.0008044386672810649, -0.0006411167399363615, -0.0004785503409299676, -0.0003172185215426642, -0.00015756701973647982, 0.0, 0.00016275479330903254, 0.0003259890748253763, 0.0004892929731254151, 0.0006522394992019593, 0.0008143920000003588]
    dvv_mwcs0_true = [-0.0008023345840610457, -0.0006390179756897862, -0.0004765897975749592, -0.00031560720282704756, -0.00015658978179569066, 0.0, 0.00016257877390881404, 0.0003256913698004309, 0.000488892319156156, 0.0006517254133130638, 0.0008137287280925279]

    @test dvv_stretch_true ≈ dvv_stretch atol=1e-16
    @test dvv_mwcs_true    ≈ dvv_mwcs atol=1e-16
    @test dvv_mwcs0_true   ≈ dvv_mwcs0 atol=1e-16
end
