using DataStructures, Dates
using SeisMonitoring: parse_inputdict

# script to make the markdown table of processing parameters

include("../../../src/Defaultproject/set_default_inputdict.jl")

InputDict_parsed = parse_inputdict(InputDict)

fo = open("../process_parameters.md", "w")

write(fo, "# Processing parameters\n\n")

write(fo, "Here is the list of processing parameters defined in the `mainparam.jl`. You can configure the parameters, which is passed to [`SeisMonitoring.run_job
`](@ref) \n\n")

write(fo, "| key | default value | type | description |\n")
write(fo, "| :--- | :--- | :--- | :--- |\n")
for (key, val) in InputDict
	txt = val[3]
	txt = replace(txt, "*" => "\\*")
	txt = replace(txt, "_" => "\\_")

	key_replaced = replace(key, "_" => "\\_")

	write(fo, "| $(key_replaced) | $(InputDict_parsed[key]) | $(val[2]) | $(txt) |\n")
end

close(fo)