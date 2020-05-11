module SeisMonitoring

# module to be used over modules
using SeisIO, SeisNoise, SeisDvv
using Dates, Printf, JLD2, Distributed, DataFrames, DataStructures
using Statistics, DSP, StatsBase
export
    seisdownload,
    seisremoveeq,
    seisxcorrelation,
    makeinput_gui,
    init_project,
    run_job






include("SeisDownload/seisdownload.jl")
include("SeisRemoveEQ/seisremoveeq.jl")
include("SeisXcorrelation/seisxcorrelation.jl")
include("SMGUI/makeinput_gui.jl")
include("Utils/init_project.jl")
include("Utils/run_job.jl")

# shared lib
include("Utils/inputdict_io.jl") #Utils, SMGUI
include("Utils/parse_inputdict.jl") #Utils, SMGUI
include("Defaultproject/set_default_inputdict.jl") # set global default InputDict

end # module
