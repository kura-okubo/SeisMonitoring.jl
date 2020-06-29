include("map_seisstack.jl")
include("seisstack_utils.jl")

"""
    seisstack(InputDict::Dict)

Compute reference (i.e. longterm) and shorttime stack of cross-correlation function.

# Major update
Instead of separating the process between stacking and dvv measurement, we perform the
SeisMeasurement duing stack in order to minimize file io cpu time.
result of dv/v, dQc-1 are attached in C.misc

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
    println("Reference starttime      = $(InputDict["reference_starttime"])")
    println("Reference endtime        = $(InputDict["reference_endtime"])")
    println("Shorttime stack starttime= $(InputDict["starttime"])")
    println("Shorttime stack endtime  = $(InputDict["endtime"])")
    println("Stacking method          = $(InputDict["stack_method"])")
    println("Compute_reference        = $(InputDict["compute_reference"])")
    println("Compute_shorttimestack   = $(InputDict["compute_shorttimestack"])")
    println("Station pairs option     = $(InputDict["stack_pairs_option"])")
    println("***************************************\n")

    # get all cc files
    cc_paths = SeisIO.ls(InputDict["stack_absolute_RawData_dir"])

    t_reference = 0.0; t_shorttime = 0.0;
    mean_assemble_cc_reference = 0.0; mean_assemble_cc_shorttime=0;mean_stack_reference=0; mean_stack_shorttime=0;
    #1. compute reference
    if InputDict["compute_reference"]
        println("-------START Reference stack--------")
        !ispath(joinpath(InputDict["fodir"], "reference")) && mkdir(joinpath(InputDict["fodir"], "reference"))

        t_reference = @elapsed bt_time_reference = pmap(x->map_seisstack(x, "reference", InputDict), cc_paths)
        # mean_assemble_cc_reference  = mean((x->x[1]).(bt_time_reference))
        # mean_stack_reference        = mean((x->x[2]).(bt_time_reference))
    end
    #2. compute shorttime stack
    if InputDict["compute_shorttimestack"]
        println("-------START Shorttime stack--------")
        !ispath(joinpath(InputDict["fodir"], "shorttime")) && mkdir(joinpath(InputDict["fodir"], "shorttime"))

        t_shorttime = @elapsed bt_time_shorttime = pmap(x->map_seisstack(x, "shorttime", InputDict), cc_paths)

        mean_assemble_cc_shorttime  = mean((x->x[1]).(bt_time_shorttime))
        mean_stack_shorttime        = mean((x->x[2]).(bt_time_shorttime))
    end

    println("seisstack has been successfully done.")

    printstyled("---Summary---\n"; color = :cyan, bold = true)
    println("Total time for reference stack = $(t_reference)[s]")
    println("Total time for shorttime stack = $(t_shorttime)[s]")
    println("mean time for reference corrdata assemble = $(mean_assemble_cc_reference)[s]")
    println("mean time for shorttime corrdata assemble = $(mean_assemble_cc_shorttime)[s]")
    println("mean time for reference stack  =$(mean_stack_reference)[s]")
    println("mean time for shorttime stack  =$(mean_stack_shorttime)[s]")

end
