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
        seisdownload::Bool=true,
        requeststation::String="",
        seisremoveeq::Bool=true,
        seisxcorrelation::Bool=true,
        seisstack::Bool=true,
        seismeasurement::Bool=true
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
    println("***************************************")



    #
    # print("Input file loading...")
    #
    #
    # println("done.")
    #

    return nothing
end
