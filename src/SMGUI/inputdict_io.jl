using Dates

function write_inputdict(f::String, InputDict::Dict)
    fo = open(f, "w")
    write(fo, "#Save at "*string(now())[1:19]*" by "*gethostname()*"\n")
    write(fo, "InputDict = Dict(\n")
    for key in keys(InputDict)
        write(fo, "\""*key*"\" => \""*InputDict[key]*"\",\n")
    end
    write(fo, ")")
    close(fo)
end
