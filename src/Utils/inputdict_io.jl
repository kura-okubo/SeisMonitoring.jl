# String in julia tuple is immutable, so tuple should be replaced with new tuple
set_values_inputdict!(InputDict, key, val) = (
    InputDict[key] = (val, InputDict[key][2], InputDict[key][3])
)

function write_inputdict(f::String, InputDict::OrderedDict)
    fo = open(f, "w")
    write(fo, "#Save at "*string(now())[1:19]*" by "*gethostname()*"\n")
    write(fo, "using Dates, DataStructures\n")
    write(fo, "InputDict = OrderedDict(\n")
    for key in keys(InputDict)
        # write(fo, "\""*key*"\" => \""*InputDict[key]*"\",\n")
        write(fo, "\""*key*"\"  => "*string(InputDict[key])*",\n")
    end
    write(fo, ")")
    close(fo)
end
