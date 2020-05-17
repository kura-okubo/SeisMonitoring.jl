using SeisMonitoring: write_inputdict
"""
    set_parameter(inputfilename::String, key::String, value::String)

Set parameter of input dictionary in the inputfile.
"""
function set_parameter(inputfilename::String, key::String, newvalue::Any)

    !ispath(inputfilename) && error("$(inputfilename) is not found in the inputfile.")
    include(inputfilename)
    !haskey(InputDict, key) && error("$(key) is not found in the inputfile.")
    #replace value
    InputDict[key] = (string(newvalue), InputDict[key][2], InputDict[key][3])
    write_inputdict(inputfilename, InputDict)

end

# set_parameter("/Users/kurama/Documents/kurama/research/testproject_1_INPUT/testparam_v1.jl", "cc_normalization", "deconvolution")
