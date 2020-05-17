module SeisMonitoring

# module to be used over modules
using SeisIO, SeisNoise, SeisDvv
using Dates, Printf, JLD2, Distributed, DataFrames, DataStructures, Distances
using Statistics, DSP, StatsBase
export

    seisdownload,
    seisremoveeq,
    seisxcorrelation,
    makeinput_gui,

    # util commands
    init_project,
    run_job,
    make_slurmbatch,
    get_parameter,
    set_parameter


include("SeisDownload/seisdownload.jl")
include("SeisRemoveEQ/seisremoveeq.jl")
include("SeisXcorrelation/seisxcorrelation.jl")
include("SeisStack/seisstack.jl")
include("SMGUI/makeinput_gui.jl")
include("Utils/init_project.jl")
include("Utils/run_job.jl")

# shared lib
include("Utils/inputdict_io.jl") #Utils, SMGUI
include("Utils/parse_inputdict.jl") #Utils, SMGUI
include("Utils/remove_nanandzerocol.jl") # used at seisstack and seismeasurement
include("Utils/slice_codawindow.jl") # used at seisstack and seismeasurement
include("Utils/get_parameter.jl") # used at seisstack and seismeasurement
include("Utils/set_parameter.jl") # used at seisstack and seismeasurement
include("Defaultproject/set_default_inputdict.jl") # set global default InputDict

# plot lib
include("SMPlot/smplot_rawdata.jl")
include("SMPlot/smplot_corrdata.jl")
include("SMPlot/smplot_stackcc.jl")

end # module
