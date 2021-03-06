include("../Defaultproject/set_default_station.jl")

"""

    init_project((
        ;
        project_name::String = "project",
        project_inputdir::String = "./",
        project_outputdir::String = "./",
        gui::String = true;
        force::Bool=false
    )

Initiate project directory where data is output.
You can output data to local machine, external HDD, scratch, etc.

# Arguments
- `project_name::String`        : project name used as directory name [default: "project"]
- `project_inputdir::String`    : absolute/relative path to make new input project directory [default: "."]
- `project_outputdir::String`   : absolute/relative path to make new output project directory [default: "."]
- `gui::Bool`                   : true if you want to use gui (you can use this to initiate case studies)
- `force::Bool=false`           : true if you want to remove existing file and init project.

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
    gui::Bool = true,
    force::Bool = false
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
        if force
            # rm(proj_input_abspath, recursive=true)
        else
            @error(
                "The project input directory already exists. Please make another project, or remove old project.",
            )
            return
        end
    end
    mkpath(proj_input_abspath)

    if ispath(proj_output_abspath)
        if force
            # rm(proj_output_abspath, recursive=true)
        else
            @error(
                "The project output directory already exists. Please make another project, or remove old project.",
            )
            return
        end
    end
    mkpath(proj_output_abspath)



    # NOTE: julia `cp` command cannot be used because of permissino issue.

    set_values_inputdict!(InputDict, "project_name", project_name)
    set_values_inputdict!(InputDict, "project_inputdir", proj_input_abspath)
    set_values_inputdict!(InputDict, "project_outputdir", proj_output_abspath)
    set_values_inputdict!(InputDict, "requeststation_file", joinpath(proj_input_abspath, "default_requeststations.jld2"))

    # make default parameter file
    write_inputdict(joinpath(proj_input_abspath, "default_param.jl"),InputDict)

    # make default request station file
    #include(joinpath(module_path, "SeisDownload/default_station.jl"))
    set_default_station(joinpath(proj_input_abspath, "default_requeststations.jld2"))
    #write_requeststation(proj_input_abspath, "request_stations.jld2", StationDataFrame)

    # make directories in OUTPUT
    mkpath(proj_output_abspath*"/seismicdata")
    mkpath(proj_output_abspath*"/cc")
    mkpath(proj_output_abspath*"/stack")
    mkpath(proj_output_abspath*"/plots")
    mkpath(proj_output_abspath*"/plots/seismicdata")
    mkpath(proj_output_abspath*"/plots/cc")
    mkpath(proj_output_abspath*"/plots/stack")
    mkpath(proj_output_abspath*"/plots/dvv")
    mkpath(proj_output_abspath*"/plots/dQc")

    println("\nProject OUTPUT $proj_output_abspath contains:")
    for (root, dirs, files) in ScanDir.walkdir(proj_output_abspath)
        for dir in dirs
            println(joinpath("/",dir)) # path to directories
        end
    end

    printstyled(
        "\n\nInput Project directory is made at:\n\n$(proj_input_abspath)\n\nOutput Project directory is made at:\n\n$(proj_output_abspath)\n\n",
        bold = true,
        color = :green,
    )

    if gui
        print("Open GUI for making input file? [y/n]: ")
        l = readline()
        if any(occursin.(["y", "Y", "yes", "Yes", "YES"], l))
            cd(proj_input_abspath)
            makeinput_gui()
        end
    end

end
