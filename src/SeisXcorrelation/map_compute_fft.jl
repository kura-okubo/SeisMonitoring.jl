"""
map_compute_fft(netstachan::String, InputDict::OrderedDict)

Compute FFT and return FFTData.

# Arguments
- `data`::Dict    : !Reprecated! Dictionary containing input data in the form: tstamp/stn => SeisChannel
- `tstamp::String,`    : Time stamp to read from JLD2 file.
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters

# Output
- `tserrorList`::Array{String,1}    : Array containing the name of failed cross-correlation timestamp/stationpairs
- `basefoname.tstamp.jld2`    : contains CorrData structure with a hierarchical structure (CC function, metadata)

# if no data is available, return empty FFTDict, which will be skipped in map_compute_cc().
"""
function map_compute_fft(netstachan::String, InputDict::OrderedDict)

    starts  = InputDict["starts_chunk"]
    ends    = InputDict["ends_chunk"]

    # NOTE: To avoid massive access to shared storage, copy files to /tmp.
    if InputDict["use_local_tmpdir"]
        local_tmp_dir = "/tmp" # default local tmp directory
        # finame = splitdir(fipath)[2]
        chunk_fi_stationdict = InputDict["chunk_fi_stationdict"]

        # return if no files exist in the rawdata path
        !haskey(chunk_fi_stationdict, netstachan) && (return (netstachan, Dict{String, FFTData}(), 0, 0))

        tmpstationlist = chunk_fi_stationdict[netstachan]
        for tmpstation in tmpstationlist
            cp(joinpath(InputDict["cc_absolute_RawData_path"], tmpstation)
                , joinpath(local_tmp_dir, tmpstation))
        end
        # fipath_tmp = joinpath(local_tmp_dir, finame)
        # cp(fipath, fipath_tmp)
        fidir = local_tmp_dir
    else
        fidir = InputDict["cc_absolute_RawData_path"]
    end

    # println("start fft processing $(netstachan).")

    t_assemble, t_fft = zeros(2)

    FFTDict = Dict{String, FFTData}()

    for tid in 1:length(starts)
        starttime, endtime = u2d.([starts[tid], ends[tid]])

        #1. assemble seisdata
        # NOTE: to avoid biased cc normalization, error if datafraction is too small (<0.95)
        # if (lowercase(InputDict["cc_normalization"]) == "coherence" ||
        #      lowercase(InputDict["cc_normalization"]) == "deconvolution") &&
        #      InputDict["data_contents_fraction"] < 0.8
        #      @warn("data_contents_fraction is better more than 0.8 with $(InputDict["cc_normalization"]) to avoid bias during spectral normalization. Please check the inputfile.")
        #  end

        # t_assemble += @elapsed S1 = assemble_seisdata_seisio(stationchannel, fi, starttime, endtime,
        #                             data_contents_fraction=InputDict["data_contents_fraction"])

        t_assemble += @elapsed S1 = assemble_seisdata(netstachan, fidir,
                                starttime, endtime, data_contents_fraction=InputDict["data_contents_fraction"])

        isnothing(S1) && continue;
        #2. applying phase shift so that the data is aligned at the sampling frequency
        phase_shift!(S1)
        #3. convert to RawData
        R1 = RawData(S1, InputDict["cc_len"], InputDict["cc_step"])
        #4. detrend, taper and band pass before computing fft and cross-correlation
        # clean_up!(R1, InputDict["freqency_band"][1], InputDict["freqency_band"][end]) # NOTE: Not applying bandpass!() here for the validation of wavelet filter
        detrend!(R1)
        demean!(R1)
        taper!(R1)

        # #5. apply one-bit normalization if true
        # InputDict["IsOnebit"] && onebit!(R1)
        #NOTE: Future work; We can make an option here to move RawData to GPU:
        #e.g. if ["GPU"]; R1 |> GPU; end
        #6. compute fft
        t_fft += @elapsed FFT1 = compute_fft(R1)
        # #7. spectral normalization
        # InputDict["cc_normalization"] == "coherence" && spectrum_coherence!(FFT1, InputDict["smoothing_windowlength"], InputDict["water_level"])
        #8. add to FFTDict
        key_FFTDict = join([netstachan, string(starttime), string(endtime)], "__") #BP.LCCB.40.SP1__2015-01-01T00:00:00__2015-01-02T00:00:00
        !isempty(FFT1) && (FFTDict[key_FFTDict] = FFT1)
    end

    if InputDict["use_local_tmpdir"]
        tmpstationlist = InputDict["chunk_fi_stationdict"][netstachan]
        for tmpstation in tmpstationlist
            rm(joinpath(local_tmp_dir, tmpstation))
        end
    end

    return (netstachan, FFTDict, t_assemble, t_fft)

end
