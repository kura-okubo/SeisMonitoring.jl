"""

    init_project(project_path::String="project")

Initiate project directory where data is output.
You can output data to local machine, external HDD, scratch, etc.

# Arguments
- `project_name`        : project name used as directory name [default: "project"]
- `project_inputdir`    : absolute/relative path to make new input project directory [default: "."]
- `project_outputdir`   : absolute/relative path to make new output project directory [default: "."]

Input and output directories can be separated for the use on local/HDD/cloud/scratch file system.

# Examples

1. run project on your local machine

```
init_project(project_name="project_test",
             project_inputdir=".",
             project_outputdir=".")
```

2. run project with external drive

```
init_project(project_name="project_external_drive",
             project_inputdir=".",
             project_outputdir="/path-to-directory-in-your-external-drive")
```
"""
function init_project(
    ;
    project_name::String = "project",
    project_inputdir::String = "./",
    project_outputdir::String = "./",
)


    proj_input_abspath = abspath(project_inputdir, project_name*"_INPUT")
    proj_output_abspath = abspath(project_outputdir, project_name*"_OUTPUT")

    if !ispath(project_inputdir)
        @error(project_inputdir * " does not exist in your system.")
        return
    elseif !ispath(project_outputdir)
        @error(project_outputdir * " does not exist in your system.")
        return
    end

    if ispath(proj_input_abspath)
        @error(
            "The project input directory already exists. Please make another project, or remove old project.",
        )
        return
    else
        mkdir(proj_input_abspath)
    end

    if ispath(proj_output_abspath)
        @error(
            "The project output directory already exists. Please make another project, or remove old project.",
        )
        return
    else
        mkdir(proj_output_abspath)
    end

    # copy default project into input path
    module_path=splitdir(pathof(SeisMonitoring))[1]

    # NOTE: julia `cp` command cannot be used because of permissino issue.
    # make default parameter file
    include(joinpath(module_path, "Defaultproject/default_param.jl"))
    include(joinpath(module_path, "SMGUI/inputdict_io.jl"))
    write_inputdict(joinpath(proj_input_abspath, "default_param.jl"),InputDict)

    # make default request station file
    #include(joinpath(module_path, "SeisDownload/default_station.jl"))
    #write_requeststation(proj_input_abspath, "request_stations.jld2", StationDataFrame)

    # make directories in OUTPUT
    mkdir(proj_output_abspath*"/seismicdata")
    mkdir(proj_output_abspath*"/cc")
    mkdir(proj_output_abspath*"/stack")
    mkdir(proj_output_abspath*"/dvv")
    mkdir(proj_output_abspath*"/dQ")
    mkdir(proj_output_abspath*"/plots")

    println("\nProject OUTPUT $proj_output_abspath contains:")
    for (root, dirs, files) in walkdir(proj_output_abspath)
        for dir in dirs
            println(joinpath("/",dir)) # path to directories
        end
    end

    printstyled(
        "\n\nInput Project directory is made at:\n\n$(proj_input_abspath)\n\nOutput Project directory is made at:\n\n$(proj_output_abspath)\n\n",
        bold = true,
        color = :green,
    )

    print("Open GUI for making input file? [y/n]: ")
    l = readline()
    if any(occursin.(["y", "Y", "yes", "Yes", "YES"], l))
        cd(proj_input_abspath)
        makeinput_gui()
    end
end
