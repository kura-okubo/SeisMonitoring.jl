module SeisMonitoring

# module to be used over modules
using SeisIO, SeisNoise, SeisDvv
using Dates, Printf, JLD2, Distributed, DataFrames, DataStructures, Distances
using Statistics, DSP, StatsBase
export

    # util commands
    init_project,
    run_job,
    makeinput_gui,
    make_slurmbatch,
    get_parameter,
    set_parameter,

    # plot command
    smplot_rawdata,
    smplot_corrdata,
    smplot_stackcc


# NOTE: Do not change the order of include.
include("SeisDownload/seisdownload.jl")
include("SeisRemoveEQ/seisremoveeq.jl")
include("SeisXcorrelation/seisxcorrelation.jl")
include("SeisStack/seisstack.jl")
include("SMGUI/makeinput_gui.jl")
include("Utils/init_project.jl")
include("Utils/run_job.jl")
include("Utils/inputdict_io.jl") #Utils, SMGUI
include("Utils/parse_inputdict.jl") #Utils, SMGUI
include("Utils/get_parameter.jl")
include("Utils/set_parameter.jl")
include("Utils/make_slurmbatch.jl")

# shared lib
include("Utils/remove_nanandzerocol.jl") # used at seisstack
include("Utils/slice_codawindow.jl") # used at seisstack
include("Utils/convert_tmpfile.jl") 
include("Defaultproject/set_default_inputdict.jl") # set global default InputDict

# plot lib
include("SMPlot/smplot_rawdata.jl")
include("SMPlot/smplot_corrdata.jl")
include("SMPlot/smplot_stackcc.jl")

end # module
