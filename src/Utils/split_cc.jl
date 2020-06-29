"""
    split_cc(x::AbstractArray, t::AbstractArray)

split cross-correlation function into positive and negative side.
"""
function split_cc(x::AbstractArray, t::AbstractArray)

    tcenter = findfirst(x -> x >= 0.0, t)
    pos_ind = tcenter:length(t)
    neg_ind = tcenter:-1:1
    t_pos   = t[pos_ind]
    t_neg   = -t[neg_ind]
    x_pos   = x[pos_ind]
    x_neg   = x[neg_ind]

    return t_pos, t_neg, x_pos, x_neg
end
