include("spectrumnormalization.jl")
include("get_cc_contents_fraction.jl")

"""
map_seisxcorrelation(station_pair::String, InputDict::OrderedDict)

Compute cross-correlation save data in jld2 file with CorrData format.

# Arguments
- `data`::Dict    : !Reprecated! Dictionary containing input data in the form: tstamp/stn => SeisChannel
- `tstamp::String,`    : Time stamp to read from JLD2 file.
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters

# Output
- `tserrorList`::Array{String,1}    : Array containing the name of failed cross-correlation timestamp/stationpairs
- `basefoname.tstamp.jld2`    : contains CorrData structure with a hierarchical structure (CC function, metadata)

#
"""
function map_seisxcorrelation(key_station_pair::String, StationPairDict::OrderedDict, InputDict::OrderedDict)

    println("start processing $(key_station_pair)")

    # open FileIO of RawData
    fi = jldopen(InputDict["cc_absolute_RawData_path"], "r")

    # open output file
    fopath = joinpath(InputDict["fodir"], key_station_pair*".jld2")
    #remove existing cc file
    ispath(fopath) && rm(fopath)
    # fo = jldopen(fopath, "w")
    # DEBUG: random save error when file size is large
    fo = jldopen(fopath, true, true, true, IOStream)

    starts  = InputDict["starts"]
    ends    = InputDict["ends"]

    #==========================#
    # 1. Compute all FFT associated with key_station_pair (e.g. "BP.CCRB-BP.EADB")
    # 2. Compute cross-correlation by the loop of station channel pairs
    #==========================#

    all_stationchannels = get_stationchanname(StationPairDict[key_station_pair])
    t_assemble, t_fft, t_xcorr = zeros(3)

    for tid in 1:length(starts)
        starttime, endtime = u2d.([starts[tid], ends[tid]])

        FFTDict = Dict()

        for stationchannel in all_stationchannels
            #1. assemble seisdata
            # NOTE: to avoid biased cc normalization, error if datafraction is too small (<0.95)
            # if (lowercase(InputDict["cc_normalization"]) == "coherence" ||
            #      lowercase(InputDict["cc_normalization"]) == "deconvolution") &&
            #      InputDict["data_contents_fraction"] < 0.8
            #      @warn("data_contents_fraction is better more than 0.8 with $(InputDict["cc_normalization"]) to avoid bias during spectral normalization. Please check the inputfile.")
            #  end

            t_assemble += @elapsed S1 = assemble_seisdata(stationchannel, fi, starttime, endtime,
                                        data_contents_fraction=InputDict["data_contents_fraction"])
            isnothing(S1) && continue;
            #2. applying phase shift so that the data is aligned at the sampling frequency
            phase_shift!(S1)
            #3. convert to RawData
            R1 = RawData(S1, InputDict["cc_len"], InputDict["cc_step"])
            #4. detrend, taper and band pass before computing fft and cross-correlation
            # clean_up!(R1, InputDict["freqency_band"][1], InputDict["freqency_band"][end]) # NOTE: Not applying bandpass!() here for the validation of wavelet filter
            detrend!(R1.x)
            taper!(R1.x,R1.fs)

            #5. apply one-bit normalization if true
            InputDict["IsOnebit"] && onebit!(R1)
            #NOTE: Future work; We can make an option here to move RawData to GPU:
            #e.g. if ["GPU"]; R1 |> GPU; end
            #6. compute fft
            t_fft += @elapsed FFT1 = compute_fft(R1)
            #7. spectral normalization
            # InputDict["cc_normalization"] == "coherence" && coherence!(FFT1, InputDict["smoothing_half_win"], InputDict["waterlevel"])
            InputDict["cc_normalization"] == "coherence" && spectrum_coherence!(FFT1, InputDict["smoothing_windowlength"], InputDict["water_level"])
            #8. add to FFTDict
            !isempty(FFT1) && (FFTDict[stationchannel] = FFT1)
        end

        # println(FFTDict)

        for channel_pair in StationPairDict[key_station_pair]
            sta1, sta2 = split(channel_pair, "-")
            # perform cross-correlation within the cc_time_unit window
            # load FFTData from dictionary
            if haskey(FFTDict, sta1)
                FFT1 = FFTDict[sta1]
            else
                println("$(sta1) at $(starttime) is missing the data. skipping cc.")
                continue
            end

            if haskey(FFTDict, sta2)
                FFT2 = FFTDict[sta2]
            else
                println("$(sta2) at $(starttime) is missing the data. skipping cc.")
                continue
            end

            # when deconvolution method id used, first station is used as source; e.g. "BP.CCRB..BP1-BP.CCRB..BP1" then BP.CCRB..BP1 is used.
            # InputDict["cc_normalization"] == "deconvolution" && deconvolution!(FFT1, InputDict["smoothing_half_win"], InputDict["waterlevel"])

            # NOTE: don't apply twice of deconvolution!(FFT1) during the loop.
            if InputDict["cc_normalization"] == "deconvolution"
                FFT1_tobecorrelated = spectrum_deconvolution(FFT1, InputDict["smoothing_windowlength"], InputDict["water_level"])
            else
                FFT1_tobecorrelated = FFT1
            end

            #8. Compute cross-correlation
            #t_xcorr += @elapsed C = compute_cc(FFT1, FFT2, InputDict["maxlag"], corr_type=InputDict["cc_method"])
            t_xcorr += @elapsed C = correlate(FFT1_tobecorrelated, FFT2, InputDict["maxlag"], corr_type=InputDict["corr_type"]) # updated version in SeisNoise.jl

            # continue if xcorr is empty
            isempty(C.corr) && continue

            # mute ccs outlier using median of maximum amplitude
            cc_medianmute!(C, InputDict["cc_medianmute_α"])

            # compute cc contents fraction
            C.misc["ccfrac_within_cc_time_unit"] = get_cc_contents_fraction(C,starttime,endtime)

            #9. Apply frequency decomposion of cross-correlation function
            C_all, freqband = compute_frequency_decomposition(C, InputDict["freqency_band"],
                                                cc_bpfilt_method=InputDict["cc_bpfilt_method"],
                                                α0=InputDict["cc_taper_α0"], αmax=InputDict["cc_taper_αmax"])

            # Prestacking to decrease the datasize if true
            # if InputDict["IsPreStack"]
            # # 1. append reference to Ctemp if IsReadReference == true for selective stack
            #     dubug_t5 += @elapsed IsReadReference && append_reference!(Ctemp, stachanpair, freqkey, ReferenceDict, InputDict)
            # # 2. perform smstack
            #     dubug_t6 += @elapsed sm_stack!(Ctemp, stackmode, InputDict) # stack with predefined stack method
            #     # println(Ctemp)
            # end

            #10. Save corr-data
            for (ic, CD) in enumerate(C_all)

                InputDict["IsPreStack"] && sm_stack!(CD, "reference", InputDict) # stack with predefined stack method

                g1 = join([string(starttime), string(endtime)], "--")  #2004-01-01T00:00:00--2004-01-02T00:00:00
                g2 = join([string(freqband[ic][1]), string(freqband[ic][2])], "-") #0.1-0.2
                groupname = joinpath(channel_pair, g1, g2)
                # println(groupname)
                # create JLD2.Group
                !haskey(fo, channel_pair) && JLD2.Group(fo, channel_pair)
                !haskey(fo[channel_pair], g1) && JLD2.Group(fo[channel_pair], g1)
                !haskey(fo, groupname) && (fo[groupname] = CD)
            end
        end
    end
    JLD2.close(fi)
    JLD2.close(fo)
    return (t_assemble, t_fft, t_xcorr)
    #DEBUG: when returning tuple on cluster, it could cause an error ProcessExitedException(57)
     #    Stacktrace:
     # [1] (::Base.var"#726#728")(::Task) at ./asyncmap.jl:178
     # [2] foreach(::Base.var"#726#728", ::Array{Any,1}) at ./abstractarray.jl:1919
     # [3] maptwice(::Function, ::Channel{Any}, ::Array{Any,1}, ::Array{String,1}) at ./asyncmap.jl:178
     # [4] wrap_n_exec_twice(::Channel{Any}, ::Array{Any,1}, ::Distributed.var"#204#207"{WorkerPool}, ::Function, ::Array{String,1}) at ./asyncmap.jl:154
     # [5] async_usemap(::Distributed.var"#188#190"{Distributed.var"#188#189#191"{WorkerPool,SeisMonitoring.var"#38#42"{OrderedCollections.OrderedDict{Any,Any},OrderedCollections.OrderedDict{String,Array{String,1}}}}}, ::Array{String,1}; ntasks::Function, batch_size::Nothing) at ./asyncmap.jl:103
     # [6] #asyncmap#710 at ./asyncmap.jl:81 [inlined]
     # [7] pmap(::Function, ::WorkerPool, ::Array{String,1}; distributed::Bool, batch_size::Int64, on_error::Nothing, retry_delays::Array{Any,1}, retry_check::Nothing) at /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.4/Distributed/src/pmap.jl:126
     # [8] pmap(::Function, ::WorkerPool, ::Array{String,1}) at /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.4/Distributed/src/pmap.jl:101
     # [9] pmap(::Function, ::Array{String,1}; kwargs::Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{,Tuple{}}}) at /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.4/Distributed/src/pmap.jl:156
     # [10] pmap at /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.4/Distributed/src/pmap.jl:156 [inlined]
     # [11] macro expansion at ./util.jl:234 [inlined]
     # [12] seisxcorrelation(::OrderedCollections.OrderedDict{String,Tuple{String,DataType,String}}) at /home1/07208/kokubo09/.julia/dev/SeisMonitoring/src/SeisXcorrelation/seisxcorrelation.jl:47
     # [13] run_job(::String; run_seisdownload::Bool, run_seisremoveeq::Bool, run_seisxcorrelation::Bool, run_seisstack::Bool) at /home1/07208/kokubo09/.julia/dev/SeisMonitoring/src/Utils/run_job.jl:119
     # [14] top-level scope at none:2
     #
    # return nothing
end
