"""

	parse_inputdict(InputDict::OrderedDict)

Parse containts of Inputdict for the use of seisdownload()
"""
function parse_inputdict(InputDict::OrderedDict)

    InputDict_parsed = OrderedDict()
    for key in keys(InputDict)

        if InputDict[key][2] == String
            # no need to parse
            InputDict_parsed[key] = InputDict[key][1]

        else
            val = parse.(InputDict[key][2], split(InputDict[key][1], ","))
            if length(val) == 1
                # this is scaler variable
                InputDict_parsed[key] = val[1]
            else
                # this is vector variables
                InputDict_parsed[key] = vec(val)
            end
        end
    end

    return InputDict_parsed
end
