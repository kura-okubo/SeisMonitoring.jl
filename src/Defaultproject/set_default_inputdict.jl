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

        "NP"                    => ("1", Int, "Number of processors you want to use for parallelization."),
        #===SeisDownload===#
        "download_time_unit"    => ("86400", Int, "[s] Unit time of data request. (e.g. request data for each 10minutes."),
        "download_margin"       => ("300", Int, "[s] Download margin to be clipped to avoid edge effect."),
        "requeststation_file"   => ("./project_OUTPUT/default_requeststations.jld2", String, "Request station dataframe saved in JLD2. See default_requeststations.jld2"),
        "savesamplefreq"        => ("20", Float64, "[Hz] Downsampling frequency. Not applied if than original sampling frequency is lower than this value."),
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

)
