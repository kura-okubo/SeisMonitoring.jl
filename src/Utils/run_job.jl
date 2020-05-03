using Distributed
using Dates

"""
    run_job(inputfile::String="";
    seisdownload::Bool=true,
    requeststation::String="",
    seisremoveeq::Bool=true,
    seisxcorrelation::Bool=true,
    seisstack::Bool=true,
    seismeasurement::Bool=true
    )

running job in the project folder.

# Arguments
- `inputfile::String`       : absolute/relative path to input file (e.g. "./project/inputfile/input.jl")

# Options
- 'seisdownload::Bool'      : run seisdownload if true [default:true]
- `requeststation::String` : absolute/relative path to request station file (e.g. "./project/inputfile/requeststation.jld2")
- 'seisremoveeq::Bool'      : run seisremoveeq if true [default:true]
- 'seisxcorrelation::Bool'  : run seisxcorrelation if true [default:true]
- 'seisstack::Bool'         : run seisstack if true [default:true]
- 'seismeasurement::Bool'   : run seismeasurement if true [default:true]
"""
function run_job(inputfile::String="";
        run_seisdownload::Bool=true,
        run_requeststation::String="",
        run_seisremoveeq::Bool=true,
        run_seisxcorrelation::Bool=true,
        run_seisstack::Bool=true,
        run_seismeasurement::Bool=true
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
    println("seisdownload       = $(seisdownload)")
    println("seisremoveeq       = $(seisdownload)")
    println("seisxcorrelation   = $(seisdownload)")
    println("seisstack          = $(seisdownload)")
    println("seismeasurement    = $(seisdownload)")
    println("***************************************\n")

    # println("Preparing parallel processing.")
    # addprocs(Inputdict["NP"])
    # println("number of procs: $(nprocs())\n")

    stall = time()

    if run_seisdownload

        st_dl=time()
        seisdownload(InputDict)
        et_dl=time()
        println("SeisDownload successfully done in $(et_dl-st_dl) seconds.\n")

    end

    if run_seisremoveeq

        st_req=time()
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

    if run_seismeasurement

        st_sm=time()
        seismeasurement(InputDict)
        et_sm=time()
        println("SeisMeasurement successfully done in $(et_sm-st_sm) seconds.\n")

    end

    etall = time()

    tall=round(etall-stall,  digits=4)
    println("*********************************************")
    println("All requested processes have been successfully done.");
    println("Computational_time[sec]:")
    print("SeisDownload, "); run_seisdownload ? println("$(et_dl-st_dl)") : println("0.0")
    print("SeisRemoveEQ, "); run_seisdownload ? println("$(et_req-st_req)") : println("0.0")
    print("SeisXcorrelation, "); run_seisdownload ? println("$(et_xc-st_xc)") : println("0.0")
    print("SeisStack, "); run_seisdownload ? println("$(et_ss-st_ss)") : println("0.0")
    print("SeisMeasurement, "); run_seisdownload ? println("$(et_sm-st_sm)") : println("0.0")
    println("Total Computational time is $(tall) seconds.");
    println("*********************************************\n")

    return nothing
end
