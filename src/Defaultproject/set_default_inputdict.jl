# Input Dictionary with Default parameters
#Format: key:variable name => (value, type, description)
#InputDict=OrderedDict{Char,Tuple}()
InputDict=OrderedDict(

        #===General parameters===#
        "project_name"          => ("project", String, "project name."),
        "project_inputdir"      => ("./project_INPUT", String, "project input directory which you initiated with init_project()."),
        "project_outputdir"     => ("./project_OUTPUT", String, "project output directory which you initiated with init_project()."),

        "starttime"             => ("2004-04-01T00:00:00", DateTime, "process start time"),
        "endtime"               => ("2004-04-03T00:00:00", DateTime, "process end time"),

        "sampling_frequency"    => ("20", Float64, "[Hz] Processing sampling frequency. Downsampling is applied when downloaded original data is higher sampling frequency."),

        "freqency_band"         => ("0.01, 0.1, 0.2, 0.5, 1.0, 2.0", Float64, "Frequency bands to be analyzed."),

        # "NP"                    => ("1", Int, "Number of processors you want to use for parallelization."),
        "MAX_MEM_USE"           => ("3.0", Float64, "[GB] Maximum memory use per core in the environment."),

        #===SeisDownload===#
        "download_time_unit"    => ("86400", Int, "[s] Unit time of data request. (e.g. request data for each 10minutes."),
        "download_margin"       => ("300", Int, "[s] Download margin to be clipped to avoid edge effect."),
        "requeststation_file"   => ("./project_OUTPUT/default_requeststations.jld2", String, "Request station dataframe saved in JLD2. See default_requeststations.jld2"),
        "IsResponseRemove"      => ("true", Bool, "True if remove instrumental response while downloading data."),
        "IsLocationBox"         => ("false", Bool, "True if using lat-lon box for request."),
        "reg"                   => ("35.7, 36.1, -120.7, -120.2", Float64, "minlat, maxlat, minlon, maxlon"),
        "Istmpfilepreserved"    => ("false", Bool, "True if you want to preserve temporal files (same size as raw data.)"),
        "IsXMLfilepreserved"    => ("false", Bool, "True if you want to preserve station xml files."),
        "numstationperrequest"  => ("1", Int, "Advanced: number of station per one HTTP request."),
        "outputformat"          => ("JLD2", String, "JLD2 or ASDF: use JLD2 if you perform the following processes with SeisMonitoring.jl"),

        #===SeisRemoveEQ===#
        "RawData_path"          => ("default", String, "\"default\" or absolute/relative path to rawdata. \"default\" links to project output directory."),
        "IsKurtosisRemoval"     => ("true", Bool, "Apply Kurtosis removal."),
        "IsSTALTARemoval"       => ("true", Bool, "Apply STA/LTA removal."),
        "IsWhitening"           => ("false", Bool, "Apply Spectral whitening."),
        "freqmin_whiten"        => ("0.1", Float64, "Minimum cutoff frequency for spectral whitening"),
        "freqmax_whiten"        => ("1.0", Float64, "Maximum cutoff frequency for spectral whitening"),
        "Append_alltraces"      => ("true", Bool, "Append kurtosis and stalta traces to SeisChannel (this increases data size)"),
        "shorttime_window"      => ("180", Float64, "Short-time window used to compute kurtosis and sta/lta"),
        "longtime_window"       => ("86400", Float64, "Long-time window used to compute sta/lta"),
        "timewindow_overlap"    => ("60", Float64, "Short-time window overlap to compute kurtosis and sta/lta"),
        "kurtosis_threshold"    => ("2.0", Float64, "Kurtosis removal threshold (The normal distribution of kurtosis is normalized to be zero.)"),
        "stalta_threshold"      => ("1.2", Float64, "STA/LTA removal threshold (For our purpose, this threshold is smaller than ordinal detection.)"),
        "stalta_absoluteclip"   => ("0.1", Float64, "[unit-of-data] clip the signal above this value (basically for instrumental error.)"),
        "fixed_tukey_margin"    => ("30", Float64, "[s] Fixed turkey margin; duration of decay outside of zero padding"),
        "IsIsolateComponents"   => ("false", Bool, "Advanced: isolating comonents at same station"),

        #===SeisXcorrelation===#
        "cc_time_unit"          => ("86400", Int64, "[s] Unit time of cross-correlation window. e.g. 60*60*24 = 86400 indicates daily-cross correlation."),
        "cc_len"                => ("3600", Int64, "[s] short-time window cross-correlation length"),
        "cc_step"               => ("1800", Int64, "[s] cross-correlation window step"),
        "maxlag"                => ("100", Float64, "[s] Maximum time lag of cross-correlation."),
        "cc_RawData_path"       => ("default", String, "\"default\" or absolute/relative path to rawdata. \"default\" links to project OUTPUT/EQRemovedData.jld2."),
        "cc_normalization"      => ("deconvolution", String, "none, coherence or deconvolution."),
        "corr_type"             => ("CC", String, "Type of correlation: `CC` (standard cross-correlation) or `PCC` (phase cross-correlation). See also doc in SeisNoise.jl"),
        "pairs_option"          => ("11, 22, 33", Array{String, 1}, "\"all\" or list of component pairs. e.g. XX, YY, ZZ"),
        "data_contents_fraction"=> ("0.8", Float64, "Advanced: discard cross-correlation if data fraction within cc_time_unit is less that this value."),
        "IsOnebit"              => ("false", Bool, "Apply One-bit normalization."),
        "smoothing_windowlength"=> ("7", Int64, "Advanced: number of points for boxcar smoothing window on coherence and deconvolution."),
        "water_level"           => ("1e-4", Float64, "Advanced: waterlevel [0.0 if not applied] on spectrum normalization with coherence and deconvolution method."),
        "cc_bpfilt_method"      => ("ButterWorth", String, "Frequency decomposition method. \"Butterworth\" or \"Wavelet\"."),
        "cc_taper_α0"           => ("0.1", Float64, "Advanced: Lowest tapering fraction for frequency adaptive tapering."),
        "cc_taper_αmax"         => ("0.25", Float64, "Advanced: Highest tapering fraction for frequency adaptive tapering."),
        "cc_medianmute_α"       => ("5.0", Float64, "Advanced: Threshold factor of median mute within cc_time_unit. NCF is removed if maximum(abs.(corr[:,i])) > cc_medianmute_α * median(maximum(abs.(corr)), dims=1)"),
        "IsPreStack"            => ("true", Bool, "Advanced: Pre-stacking corrdata within each cc_time_unit when assembling the corrdata for the sake of saving memory use."),

        #===SeisStack===#
        "stack_RawData_dir"     => ("default", String, "\"default\" or absolute/relative path to cc directory. \"default\" links to project OUTPUT/cc."),
        "stack_method"          => ("linear", String, "stacking method: linear, selective, robust, pws, robustpws are available"),
        "compute_reference"     => ("true", Bool, "true if compute reference stack for longterm stack."),
        "compute_shorttimestack"=> ("true", Bool, "true if compute shorttime stack for continuous monitoring."),
        "stack_pairs_option"    => ("11, 22, 33", Array{String, 1}, "\"all\" or list of component pairs. e.g. XX, YY, ZZ"),
        "averagestack_factor"   => ("1", Int, "Integer factor of cc_time_unit for stacking duration. e.g. cc_time_unit = 1day and averagestack_factor=30 provides 30days moving window average."),
        "averagestack_step"     => ("1", Int, "Step of averagestack window."),
        "min_cc_datafraction"   => ("0.5", Float64, "Advanced: discard cross-correlation if data fraction within stacking period is less that this value."),
        "reference_starttime"   => ("2004-04-01T00:00:00", DateTime, "reference start time"),
        "reference_endtime"     => ("2004-04-02T00:00:00", DateTime, "reference end time"),
        "dist_threshold"        => ("1.0", Float64, "Threshold of distance used for selective stacking."),
        "distance_type"         => ("CorrDist", String, "Advanced: Distance type used in selective stacking. See https://github.com/JuliaStats/Distances.jl for available types."),
        "IsZeropadBeforeStack"  => ("false", Bool, "Zero padding outside of coda window using tukey window before stacking."),
        "background_vel"        => ("2000.0", Float64, "[m/s] Approximation of background wave velocity, just used for coda slicing."),
        "coda_Qinv"             => ("0.01", Float64, "Approximation of inverse coda_Q, Qc^{-1}, just used for coda slicing."),
        "min_ballistic_twin"    => ("5.0", Float64, "[s] Explicit ballistic time window to remove coherence around zero timelag. This is aimed to remove it mainly for auto-correlation."),
        "max_coda_length"       => ("60.0", Float64, "[s] Maximum coda window length."),
        "slice_minthreshold"    => ("0.1", Float64, "Advanced: Threshold for attenuation decay."),
        "IsAlternateRefChannel" => ("true", Bool, "Advanced: Allow for using alternative station channel for reference. (e.g. BP.LCCB..BP1-BP.MMNB..BP1 is used as reference for BP.LCCB..SP1-BP.MMNB..SP1)"),

        #SeisMeasurement
        "measurement_method"    => ("dualstretching", String, "Stretching method for measuring dv/v and dQ^{-1}. \"stretching\",\"mwcs\",\"wcc\",\"dtw\",\"dualstretching\" "),
        "dvv_stretching_range"  => ("0.02", Float64, "Advanced: dvv stretching trial range for dvv (+- abs(dvv_stretching_range))."),
        "dvv_stretching_Ntrial" => ("101", Int, "Advanced: dvv stretching trial number for dvv."),
        "dQc_stretching_range"  => ("0.004", Float64, "Advanced: dual stretching trial range for dQc (+- abs(dvv_stretching_range))."),
        "dQc_stretching_Ntrial" => ("51", Int, "Advanced: dual stretching trial number for dQc."),
        "dAA_stretching_range"  => ("0.05", Float64, "Advanced: dual stretching trial range for dAA (+- abs(dvv_stretching_range))."),
        "dAA_stretching_Ntrial" => ("51", Int, "Advanced: dual stretching trial number for dAA."),
        "smoothing_window_len"  => ("15.0", Float64, "Advanced: [s] smoothing time window length to measure coda Q change."),
        "stretch_distmethod"    => ("euclidean", String, "Distance method used for stretching method. \"euclidean\" or \"mahalanobis\""),
        "stretch_debugplot"     => ("true", Bool, "Advanced: plot debug figures for stretching method (this option increases stacking cpu time)"),

)
