include("map_seisstack.jl")
include("seisstack_utils.jl")
include("assemble_corrdata.jl")

"""
    seisstack(InputDict::Dict)

Compute reference (i.e. longterm) and shorttime stack of cross-correlation function.

# Arguments
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters
"""
function seisstack(InputDict_origin::OrderedDict)

    InputDict = parse_inputdict(InputDict_origin)

    project_outputdir = abspath(InputDict["project_outputdir"])
    InputDict["fodir"] = joinpath(project_outputdir, "stack")

    if lowercase(InputDict["stack_RawData_dir"]) == "default"
        InputDict["stack_absolute_RawData_dir"] = joinpath(project_outputdir, "cc")
    else
        InputDict["stack_absolute_RawData_dir"] = abspath(InputDict["stack_RawData_dir"])
    end

    println("***************************************")
    println("Stacking starttime= $(InputDict["starttime"])")
    println("Stacking endtime = $(InputDict["endtime"])")
    println("Stacking method = $(InputDict["stack_method"])")
    println("Compute_reference = $(InputDict["compute_reference"])")
    println("Compute_shorttimestack = $(InputDict["compute_shorttimestack"])")
    println("Station pairs option = $(InputDict["stack_pairs_option"])")
    println("***************************************\n")

    # get all cc files
    cc_paths = SeisIO.ls(InputDict["stack_absolute_RawData_dir"])

    t_reference = 0.0
    t_shorttime = 0.0
    bt_time_reference = 0.0
    bt_time_shorttime = 0.0
    #1. compute reference
    if InputDict["compute_reference"]
        println("-------START Reference stack--------")
        mkdir(joinpath(InputDict["fodir"], "reference"))

        t_reference = @elapsed bt_time_reference = pmap(x->map_seisstack(x, "reference", InputDict), cc_paths)
    end
    #2. compute shorttime stack
    if InputDict["compute_shorttimestack"]
        println("-------START Shorttime stack--------")
        mkdir(joinpath(InputDict["fodir"], "shorttime"))

        t_shorttime = @elapsed bt_time_shorttime = pmap(x->map_seisstack(x, "shorttime", InputDict), cc_paths)
    end

    println("seisstack has been successfully done.")

    mean_assemble_cc_reference   = mean((x->x[1]).(bt_time_reference))
    mean_assemble_cc_shorttime   = mean((x->x[1]).(bt_time_shorttime))
    mean_stack_reference        = mean((x->x[2]).(bt_time_reference))
    mean_stack_shorttime        = mean((x->x[2]).(bt_time_shorttime))

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("Total time for reference stack = $(t_reference)[s]")
    println("Total time for shorttime stack = $(t_shorttime)[s]")
    println("mean time for reference corrdata assemble = $(mean_assemble_cputime)[s]")
    println("mean time for shorttime corrdata assemble = $(mean_fft_cputime)[s]")
    println("mean time for reference stack  =$(mean_stack_reference)[s]")
    println("mean time for shorttime stack  =$(mean_stack_shorttime)[s]")

end
