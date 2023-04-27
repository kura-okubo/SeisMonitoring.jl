#including stacking method
include("selectivestack.jl")

function sm_stack!(C::CorrData, stackmode::String, InputDict::OrderedDict)

    stack_method = InputDict["stack_method"]
    t_beforestack = C.t

    if lowercase(stack_method) == "linear"
        SeisNoise.stack!(C, allstack=true, stacktype=mean)

    elseif lowercase(stack_method) == "selective"
        selectivestack!(C, stackmode, dist_threshold=InputDict["dist_threshold"],
                        distance_type=InputDict["distance_type"])

    elseif lowercase(stack_method) == "robust"
        SeisNoise.robuststack!(C)

    elseif lowercase(stack_method) == "pws"
        SeisNoise.pws!(C)

    elseif lowercase(stack_method) == "robustpws"
        SeisNoise.robustpws!(C)

    else
        error("stack_method: $(stack_method) is not available.")
    end

    C.misc["t_stacked"] = t_beforestack
    C.misc["stack_method"] = stack_method
end
