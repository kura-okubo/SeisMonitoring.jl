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

        "freqency_band"         => ("0.01, 0.1, 0.1, 0.2", Float64, "Frequency bands to be analyzed.")

        "NP"                    => ("1", Int, "Number of processors you want to use for parallelization."),
        #===SeisDownload===#
        "download_time_unit"    => ("86400", Int, "[s] Unit time of data request. (e.g. request data for each 10minutes."),
        "download_margin"       => ("300", Int, "[s] Download margin to be clipped to avoid edge effect."),
        "requeststation_file"   => ("./project_OUTPUT/default_requeststations.jld2", String, "Request station dataframe saved in JLD2. See default_requeststations.jld2"),
        "IsResponseRemove"      => ("true", Bool, "True if remove instrumental response while downloading data."),
        "IsLocationBox"         => ("false", Bool, "True if using lat-lon box for request."),
        "reg"         => ("35.7, 36.1, -120.7, -120.2", Float64, "minlat, maxlat, minlon, maxlon"),
        "Istmpfilepreserved"    => ("false", Bool, "True if you want to preserve temporal files (same size as raw data.)"),
        "IsXMLfilepreserved"    => ("false", Bool, "True if you want to preserve station xml files."),
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
        "cc_method"             => ("cross-correlation", String, "cross-correlation, coherence or deconvolution."),
        "pairs_option"          => ("11, 22, 33", Array{String, 1}, "\"all\" or list of component pairs. e.g. XX, YY, ZZ"),
        "data_contents_fraction"=> ("0.8", Float64, "Advanced: discard cross-correlation if data fraction within cc_time_unit is less that this value."),
        "IsOnebit"              => ("false", Bool, "Apply One-bit normalization."),

)
