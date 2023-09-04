module SeisMonitoring

# module to be used over modules
using SeisIO, SeisNoise, SeisDvv
using Dates, Printf, JLD2, Distributed, DataFrames, DataStructures, Distances
using Statistics, DSP, StatsBase, ScanDir
export

    # util commands
    init_project,
    run_job,
    # makeinput_gui,
    make_slurmbatch,
    make_requeststation_fromIRISgmap,
    get_parameter,
    set_parameter,

    # plot command
    smplot_rawdata,
    smplot_corrdata,
    smplot_stackcc,
    smplot_noiseavailability,
    smplot_pdfpsd,
    smplot_pdfdvv,

    #stats command
    smstats_dataavailability,
    smstats_read_computedvvdqq,
    smstats_read_mwcs,
    smstats_read_stretching

# NOTE: Do not change the order of include.
include("SeisDownload/seisdownload.jl")
include("SeisRemoveEQ/seisremoveeq.jl")
include("SeisXcorrelation/seisxcorrelation.jl")
include("SeisStack/seisstack.jl")
include("SeisStack/sm_stack.jl") # used both map_seisstack and assemble_corrdata for prestacking
# NOTE: do not include GUI if julia is in parallel.
# include("SMGUI/makeinput_gui.jl")
include("Utils/init_project.jl")
include("Utils/run_job.jl")
include("Utils/inputdict_io.jl") #Utils, SMGUI
include("Utils/parse_inputdict.jl") #Utils, SMGUI
include("Utils/get_parameter.jl")
include("Utils/set_parameter.jl")
include("Utils/make_slurmbatch.jl")
include("Utils/make_requeststation.jl")

# shared lib
include("Utils/remove_nanandzerocol.jl") # used at seisstack
# include("Utils/slice_codawindow.jl") # used at seisstack
# include("Utils/energybased_slice_codawindow.jl") # used at seisstack
# include("Utils/log10_slice_codawindow.jl") # used at seisstack
# include("Utils/mwcc_slice_codawindow.jl") # used at seisstack
include("Utils/const_slice_codawindow.jl") # used at seisstack
include("Utils/convert_tmpfile.jl") # used at seisdownload and seisremoveeq
include("Utils/get_noisedatafraction.jl") # used at seisdownload, seisremoveeq and seisxcorrelation
include("Utils/split_cc.jl") # used at slice_coda_window, seismeasurement
include("Utils/get_chanpairtype.jl") # used at seisxcorrelation, seisstack

include("Defaultproject/set_default_inputdict.jl") # set global default InputDict

# plot lib
include("SMPlot/smplot_rawdata.jl")
include("SMPlot/smplot_corrdata.jl")
include("SMPlot/smplot_stackcc.jl")
include("SMPlot/smplot_noiseavailability.jl")
include("SMPlot/smplot_ppsd.jl")
include("SMPlot/smplot_pdfdvv.jl")
include("SMStats/smstats_dataavailability.jl")
include("SMStats/smstats_read_computedvvdqq.jl")
include("SMStats/smstats_read_stretching.jl")
include("SMStats/smstats_read_mwcs.jl")

end # module
