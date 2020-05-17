using SeisMonitoring: parse_inputdict

"""
    make_slurmbatch(inputfile)

Submit jobfile to HPC cluster using slurm batch file and execute run_job().

# Argument

- `inputfile::String`       : **Absolute** path to input file (e.g. "./project/inputfile/input.jl")

# Options
- 'run_seisdownload::Bool'      : run seisdownload if true [default:true]
- 'run_seisremoveeq::Bool'      : run seisremoveeq if true [default:true]
- 'run_seisxcorrelation::Bool'  : run seisxcorrelation if true [default:true]
- 'run_seisstack::Bool'         : run seisstack if true [default:true]
"""
function make_slurmbatch(inputfile::String="",
        partition::String="normal",
        totalnodes::Int=1,
        totaltasks::Int=32,
        runtime::String="01:00:00", # write in slurm format
        absolute_juliapath::String="/Applications/Julia-1.4.app/Contents/Resources/julia/bin/julia"; # please find this by type `which julia` in bash shell
        run_seisdownload::Bool=true,
        run_seisremoveeq::Bool=true,
        run_seisxcorrelation::Bool=true,
        run_seisstack::Bool=true,
        run_seismeasurement::Bool=true
        )


        ispath(inputfile) ? include(inputfile) : error("$(inputfile) not found.")
        InputDict_parsed = parse_inputdict(InputDict)
        project_name = InputDict_parsed["project_name"]
        project_inputdir  = InputDict_parsed["project_inputdir"]
        project_outputdir = InputDict_parsed["project_outputdir"]

        batchfile =  joinpath(project_inputdir, "slurm_run_seismonitoring.slurm")
        out = joinpath(project_outputdir, "out_$(project_name)_%j.txt")
        err = joinpath(project_outputdir, "err_$(project_name)_%j.txt")

        fo = open(batchfile, "w")
        write(fo, "#!/bin/bash\n#----------------------------------------------------\n")
        write(fo, "#SBATCH -p $(partition)\n")
        write(fo, "#SBATCH -j $(project_name)\n")
        write(fo, "#SBATCH -N $(totalnodes)\n")
        write(fo, "#SBATCH -n $(totaltasks)\n")
        write(fo, "#SBATCH -t $(runtime)\n")
        write(fo, "#SBATCH -o $(out)\n")
        write(fo, "#SBATCH -e $(err)\n")
        write(fo, "#----------------------------------------------------\n\n")

        write(fo, "date\n")
        # write(fo, "module purge\n")
        write(fo, "xvfb-run $(absolute_juliapath) -p $(totaltasks) -e 'using SeisMonitoring;
        run_job(\"$(inputfile)\", run_seisdownload=$(run_seisdownload),
        run_seisremoveeq=$(run_seisremoveeq),
        run_seisxcorrelation=$(run_seisxcorrelation),
        run_seisstack=$(run_seisstack))'\n")
        write(fo, "date\n")
        write(fo, "echo \"job $(project_name) is successfully done.\"")
        close(fo)

        @info("$(batchfile) is successfully made.
Please exit julia, and type shell command:
`sbatch $(batchfile)`
in your terminal.")
end
