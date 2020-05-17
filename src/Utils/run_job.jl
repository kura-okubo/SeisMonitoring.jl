using Distributed

include("printlogos.jl")

"""
    run_job(inputfile::String="";
    seisdownload::Bool=true,
    seisremoveeq::Bool=true,
    seisxcorrelation::Bool=true,
    seisstack::Bool=true,
    )

running job in the project folder.

# Arguments
- `inputfile::String`       : absolute/relative path to input file (e.g. "./project/inputfile/input.jl")

# Options
- 'seisdownload::Bool'      : run seisdownload if true [default:true]
- 'seisremoveeq::Bool'      : run seisremoveeq if true [default:true]
- 'seisxcorrelation::Bool'  : run seisxcorrelation if true [default:true]
- 'seisstack::Bool'         : run seisstack if true [default:true]
"""
function run_job(inputfile::String="";
        run_seisdownload::Bool=true,
        run_seisremoveeq::Bool=true,
        run_seisxcorrelation::Bool=true,
        run_seisstack::Bool=true,
        )

    if isempty(inputfile)
        @error("Please enter the input parameter file name. (e.g. ./default_param.jl)")
        return
    end

    inputfile_abspath=abspath(inputfile)
    if !ispath(inputfile_abspath)
        @error("$(inputfile_abspath) is not found.")
        return
    end

    initlogo()

    printstyled(
        "\ninput file: $(inputfile_abspath)\n\n",
        bold = true,
        color = :cyan,
    )

    println("***************************************")
    println("seisdownload       = $(run_seisdownload)")
    println("seisremoveeq       = $(run_seisremoveeq)")
    println("seisxcorrelation   = $(run_seisxcorrelation)")
    println("seisstack          = $(run_seisstack)")
    println("***************************************\n")

    include(inputfile)

    println("Preparing parallel processing.")
    #=== NOTE: Julia `nprocs` is not the number of processors to be used.
    Julia nprocs is the number of processes to be parallelized using available
    processors in the system.
    However, for the sake of simplicity, we indicate np as number of processors
    following traditional use and addprocs according to the np.
    So you can assign the large number of NP; it is just parallelized by the number
    of available cores in your system.
    ===#
    procs_tobeadded = parse(Int, InputDict["NP"][1]) - nprocs()
    if procs_tobeadded < 0
        rmprocs(abs(procs_tobeadded))
    elseif procs_tobeadded >= 1
        addprocs(procs_tobeadded)
    end
    println("NP               : $(nprocs())")
    if nprocs()-1 >= 1; println("Number of workers: $(nprocs()-1)\n") end
    eval(macroexpand(SeisMonitoring, quote @everywhere using SeisMonitoring end))

    stall = time()

    if run_seisdownload

        st_dl=time()
        printstyled(
            "\nStart running SeisDownload\n\n",
            bold = true,
            color = :cyan,
        )
        seisdownload(InputDict)
        et_dl=time()
        println("SeisDownload successfully done in $(et_dl-st_dl) seconds.\n")

    end

    if run_seisremoveeq

        st_req=time()
        printstyled(
            "\nStart running SeisRemoveEQ\n\n",
            bold = true,
            color = :cyan,
        )
        seisremoveeq(InputDict)
        et_req=time()
        println("SeisRemoveEQ successfully done in $(et_req-st_req) seconds.\n")

    end

    if run_seisxcorrelation

        st_xc=time()
        seisxcorrelation(InputDict)
        et_xc=time()
        println("SeisXcorrelation successfully done in $(et_xc-st_xc) seconds.\n")

    end

    if run_seisstack

        st_ss=time()
        seisstack(InputDict)
        et_ss=time()
        println("SeisStack successfully done in $(et_ss-st_ss) seconds.\n")

    end

    etall = time()

    tall=round(etall-stall,  digits=4)
    println("*********************************************")
    println("All requested processes have been successfully done.");
    println("Computational_time[sec]:")
    print("SeisDownload, "); run_seisdownload ? println("$(et_dl-st_dl)") : println("0.0")
    print("SeisRemoveEQ, "); run_seisremoveeq ? println("$(et_req-st_req)") : println("0.0")
    print("SeisXcorrelation, "); run_seisxcorrelation ? println("$(et_xc-st_xc)") : println("0.0")
    print("SeisStack, "); run_seisstack ? println("$(et_ss-st_ss)") : println("0.0")
    println("Total Computational time is $(tall) seconds.");
    println("*********************************************\n")

    return nothing
end
