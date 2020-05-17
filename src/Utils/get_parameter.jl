using Dates, DataStructures

"""
    get_parameter(inputfilename::String, parametername::String)

get parameter from input file.
"""
function get_parameter(inputfilename::String, key::String)
    !ispath(inputfilename) && error("$(inputfilename) is not found.")
    include(inputfilename)
    !haskey(InputDict, key) && error("$(key) is not found.")
    if InputDict[key][2] == String
        # no need to parse
        return InputDict[key][1]

    elseif InputDict[key][2] == Array{String, 1}
        # array of string e.g. "XX, YY, ZZ"
        return lstrip.(split(InputDict[key][1], ","))

    else
        val = parse.(InputDict[key][2], split(InputDict[key][1], ","))
        if length(val) == 1
            # this is scaler variable
            return val[1]
        else
            # this is vector variables
            return vec(val)
        end
    end
end

# finame = "/Users/kurama/Documents/kurama/research/testproject_1_INPUT/testparam_v1.jl"
#
# get_parameter(finame, "project_outputdir")
# get_parameter(finame, "cc_time_unit")
# get_parameter(finame, "stack_pairs_option")
