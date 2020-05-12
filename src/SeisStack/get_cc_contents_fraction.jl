"""
    get_cc_contents_fraction()

return cross-correlation coverage between starttime and endtime

# Algorithm
Iterating to merge adjusent time windows. When all time windows are merged,
compute total duration of data and its fraction to the start and end time.

# NOTE
This can be used **only before stacking** as it assumes `[t-cc_len/2, t+cc_len/2]`.
"""
function get_cc_contents_fraction(
    C::CorrData,
    starttime::DateTime,
    endtime::DateTime
)
    # initiate current os and oe; ot_cur = [os_cur, oe_cur]
    ot_cur = map(i-> [C.t[i] - C.cc_len/2, C.t[i] + C.cc_len/2], 1:length(C.t))
    ot_merged = Array{Array{Float64,1}, 1}(undef, 0)

    # if C.t is only one component, skip iteration process
    if isone(length(ot_cur))
        ot_merged = ot_cur
    else
        while true
            s1, e1 = ot_cur[1]
            s2, e2 = ot_cur[2]
            ot = t_coverage(s1, e1, s2, e2)
            if iszero(ot)
                # this time window is isolated. add to ot_merged and popfirst!()
                push!(ot_merged, [s1, e1])
                popfirst!(ot_cur)
            else
                # the first and second row of t is merged
                popfirst!(ot_cur)
                popfirst!(ot_cur)
                pushfirst!(ot_cur, ot)
            end

            if isone(length(ot_cur))
                # this is the last chunk of merged time window
                push!(ot_merged, ot_cur[1])
                break;
            end
        end
    end
    os = (x->x[1]).(ot_merged)
    oe = (x->x[2]).(ot_merged)
    # compute duration of cc time windows
    cc_contents_total = sum(oe .- os)
    target_window_length = (endtime - starttime).value / 1e3 # millisecond -> second
    return cc_contents_total / target_window_length
end

"""
    t_coverage(s1::Int, e1::Int, s2::Int, e2::Int)
return time coverage between two time window
# Argument
-`s1,e1,s2,e2`: Start and end time in unix time [s]
"""
function t_coverage(s1::AbstractFloat, e1::AbstractFloat, s2::AbstractFloat, e2::AbstractFloat)
    if s2>e1 || s1>e2
        # no overlap
        return 0
    else
        os = min(s1, s2)
        oe = max(e1, e2)
        return [os, oe]
    end
end
#
# c1 = jldopen("corr_test.jld2", "r")
# C = c1["test"]
# close(c1)
# # manipulate time
# starttime = DateTime("2004-04-01T00:00:00")
# endtime = DateTime("2004-04-02T00:00:00")
#
# #tinds = vcat(collect(1:2:10), collect(20:1:46))
# # tinds = vcat(collect(1:5:20))
# tinds = 1:1
# C.t = C.t[tinds]
# C.corr = C.corr[:, tinds]
#
# cc_data_fraction = get_cc_contents_fraction(C,starttime,endtime)
