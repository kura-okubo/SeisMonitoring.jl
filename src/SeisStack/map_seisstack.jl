include("assemble_corrdata.jl")
include("seismeasurement.jl")
#including stacking method
include("selectivestack.jl")

"""
    map_seisstack(fipath, stackmode::String, InputDict::OrderedDict)

Stack cross-correlation function
"""
function map_seisstack(fipath, stackmode::String, InputDict::OrderedDict)

    println("start processing $(splitdir(fipath)[2][1:end-5])")

    # evaluate if we need to read and append reference
    IsReadReference = (stackmode=="shorttime")

    # open FileIO of cross-correlation
    fi = jldopen(fipath, "r")

    # open output file
    if stackmode=="reference"
        foname = "reference_"*splitdir(fipath)[2]
        fopath = joinpath(InputDict["fodir"], "reference", foname)
        # single start-end time for reference
        starts, ends = [[InputDict["reference_starttime"]], [InputDict["reference_endtime"]]]
    elseif stackmode=="shorttime"
        foname = "shorttime_"*splitdir(fipath)[2]
        fopath = joinpath(InputDict["fodir"], "shorttime", foname)
        # multiple shorttime stacking windows
        starts, ends = get_shorttime_window(InputDict["starttime"],InputDict["endtime"], InputDict["cc_time_unit"],
                                                    InputDict["averagestack_factor"], InputDict["averagestack_step"])
    end

    ispath(fopath) && rm(fopath)
    fo = jldopen(fopath, "w")

    #==========================#
    # 1. Assemble corrdata between starttime and endtime
    # 2. Compute stacking
    #==========================#

    # Reading reference traces if stack_method needs reference
    if IsReadReference
        # load reference
        refname = "reference_"*splitdir(fipath)[2] #reference_BP.CCRB-BP.EADB.jld2
        fi_refpath =  joinpath(InputDict["fodir"], "reference", refname)
        if ispath(fi_refpath)
            fi_ref = jldopen(fi_refpath, "r")
            @show ReferenceDict = get_reference(fi_ref)
            close(fi_ref)
        else
            # No reference found for this station pairs
            ReferenceDict = Dict()
        end
    end

    t_assemblecc = 0; t_stack = 0; t_seismeasurement = 0;

    for stachanpair = collect(keys(fi))

        # stachanpair: BP.CCRB..BP1-BP.EADB..BP1
        sta1, sta2 = split(stachanpair, "-")
        comp = sta1[end]*sta2[end]

        # skip if component pair is not in the list
        (comp ∉ InputDict["stack_pairs_option"] && InputDict["stack_pairs_option"] != "all") && continue;

        CorrData_Buffer = Dict()

        for tid in 1:length(starts)
            starttime, endtime = [starts[tid], ends[tid]]
            centraltime = u2d((d2u(starttime) + d2u(endtime)) /2) #central time between starttime and endtime

            println("start processing $(stachanpair) at $(string(starttime))-$(string(endtime))")

            # assemble corrdata
            t_assemblecc += @elapsed C_all, CorrData_Buffer = assemble_corrdata(fi,stachanpair,starttime,endtime,InputDict["freqency_band"],
                                    CorrData_Buffer=CorrData_Buffer,
                                    min_cc_datafraction = InputDict["min_cc_datafraction"],
                                    MAX_MEM_USE=InputDict["MAX_MEM_USE"])

            # stack with respect to frequency band
            for freqkey in collect(keys(C_all))
                C = C_all[freqkey]

                (isempty(C.corr) || isempty(C.t)) && continue  # this does not have cc trace within the time window.
                remove_nanandzerocol!(C)  # remove column which has NaN or all zero
                (isempty(C.corr) || isempty(C.t)) && continue  # this does not have cc trace within the time window.

                # slice coda window and zero padding before stack if true
                coda_window, timelag, fillbox = slice_codawindow!(C,
                                        InputDict["background_vel"],
                                        InputDict["coda_Qinv"],
                                        InputDict["min_ballistic_twin"],
                                        InputDict["max_coda_length"],
                                        attenuation_minthreshold=InputDict["slice_minthreshold"],
                                        zeropad=InputDict["IsZeropadBeforeStack"])
                # append coda_window, timelag and fillbox for plotting
                C.misc["coda_window"] = coda_window
                C.misc["timelag"] = timelag
                C.misc["fillbox"] = fillbox

                # append reference curve if needed
                IsReadReference && append_reference!(C, stachanpair, freqkey, ReferenceDict, InputDict)

                t_stack += @elapsed sm_stack!(C, stackmode, InputDict) # stack with predefined stack method

                (isempty(C.corr) || isempty(C.t)) && continue  # this does not have cc trace within the time window.
                remove_nanandzerocol!(C)  # remove column which has NaN or all zero
                (isempty(C.corr) || isempty(C.t)) && continue  # this does not have cc trace within the time window.

                # append metadata
                C.misc["stack_starttime"] = starttime
                C.misc["stack_endtime"] = endtime
                C.misc["stack_centraltime"] = centraltime

                # compute dv/v and dQinv
                stackmode=="shorttime" && (t_seismeasurement += @elapsed seismeasurement!(C, InputDict))

                # save data into jld2
                g1 = join([string(starttime), string(endtime)], "--")  #2004-01-01T00:00:00--2004-01-02T00:00:00
                g2 = freqkey #0.1-0.2
                groupname = joinpath(stachanpair, g1, g2)
                # println(groupname)
                # create JLD2.Group
                !haskey(fo, stachanpair) && JLD2.Group(fo, stachanpair)
                !haskey(fo[stachanpair], g1) && JLD2.Group(fo[stachanpair], g1)
                !haskey(fo, groupname) && (fo[groupname] = C)
            end

            println("debug: current t_seismeasurement = $(t_seismeasurement)[s]")
        end
    end

    close(fi)
    close(fo)
    return (t_assemblecc, t_stack, t_seismeasurement)
end

function sm_stack!(C::CorrData, stackmode::String, InputDict::OrderedDict)

    stack_method = InputDict["stack_method"]

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

    C.misc["stack_method"] = stack_method
end
