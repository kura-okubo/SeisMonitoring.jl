# Input Dictionary with Default parameters
#Format: key:variable name => (value, type, description)
#InputDict=OrderedDict{Char,Tuple}()
InputDict=OrderedDict(

        #===General parameters===#
        "project_name"          => ("project", String, "project name."),
        "project_inputdir"      => ("./project_INPUT", String, "project input directory which you initiated with init_project()."),
        "project_outputdir"     => ("./project_OUTPUT", String, "project output directory which you initiated with init_project()."),

        "starttime"             => ("2020-04-01T00:00:00", DateTime, "process start time"),
        "endtime"               => ("2020-04-11T00:00:00", DateTime, "process end time"),

        "NP"                    => ("1", Int, "Number of processors you want to use for parallelization."),
        #===SeisDownload===#
        "download_time_unit"    => ("600", Int, "[s] Unit time of data request. (e.g. request data for each 10minutes."),
        "download_margin"       => ("300", Int, "[s] Download margin to be clipped to avoid edge effect."),
        "requeststation_file"   => ("./project_OUTPUT/default_requeststations.jld2", String, "Request station dataframe saved in JLD2. See default_requeststations.jld2"),
        "savesamplefreq"        => ("20", Float64, "[Hz] Downsampling frequency. Not applied if than original sampling frequency is lower than this value."),
        "IsResponseRemove"      => ("true", Bool, "True if remove instrumental response while downloading data."),
        "IsLocationBox"         => ("false", Bool, "True if using lat-lon box for request."),
        "reg"         => ("35.7, 36.1, -120.7, -120.2", Float64, "minlat, maxlat, minlon, maxlon"),
        "Istmpfilepreserved"    => ("false", Bool, "True if you want to preserve temporal files (same size as raw data.)"),
        "IsXMLfilepreserved"    => ("false", Bool, "True if you want to preserve station xml files."),
        "outputformat"          => ("JLD2", String, "JLD2 or ASDF: use JLD2 if you perform the following processes with SeisMonitoring.jl"),
)
