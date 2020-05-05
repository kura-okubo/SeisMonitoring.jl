include("seisremoveeq_utils.jl")
include("map_removeEQ.jl")

"""
	seisremoveeq(InputDict::OrderedDict)

remove transient signals using STA/LTA and Kurtosis
"""
function seisremoveeq(InputDict_origin::OrderedDict)

	InputDict = parse_inputdict(InputDict_origin)

	project_outputdir		= abspath(InputDict["project_outputdir"])
	InputDict["fodir"] 		= joinpath(project_outputdir, "seismicdata")
	tmpdir 					= joinpath(project_outputdir, "seismicdata", "seisremoveeq_tmp")
	InputDict["tmpdir_rem"] = tmpdir

	if ispath(tmpdir); rm(tmpdir, recursive=true); end
	mkdir(tmpdir)

	# parallelize with keys in Rawdata.jld2 i.e. stations
	if RawData_path = "default"
		RawData_path = joinpath(InputDict["fodir"], "RawData.jld2")
	else
		RawData_path = InputDict["RawData_path"]
	end

	t = jldopen(RawData_path, "r")

	t_removeeq = @elapsed pmap(x -> map_removeEQ(x, t, InputDict), keys(t))
	t_convert = @elapsed convert_tmpfile(InputDict)

	printstyled("---Summary---\n"; color=:cyan, bold=true)
	println("time to remove EQ    =$(t_removeeq)[s]")
	println("time to convert      =$(t_convert)[s]")

end
