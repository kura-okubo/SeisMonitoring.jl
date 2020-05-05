__precompile__()
module Get_kurtosis

export get_kurtosis

using Statistics, SeisIO

"""
    get_kurtosis(data::SeisChannel,kurtsis_tw_sparse::Float64; timewinlength::Float64=60)

    compute kurtosis at each timewindow

# Input:
    - `data::SeisChannel`    : SeisData from SeisIO
    - `kurtosis_tw_sparse::Float64` : time length of span for kurtosis time window
    - `timewinlength::Float64`  : time window to calculate kurtosis
    kurtosis evaluation following Baillard et al.(2013)
"""
function get_kurtosis(data::SeisChannel, timewinlength::Float64=180.0, kurtosis_tw_sparse::Float64=60.0)

    #convert window lengths from seconds to samples
    TimeWin = trunc(Int, timewinlength * data.fs)
    SparseWin = trunc(Int, kurtosis_tw_sparse * data.fs)
    data.misc["kurtosis"] = fast_kurtosis_series(data.x, TimeWin, SparseWin)

    return data

end


"""
    fast_kurtosis_series(v::RealArray, TimeWin::Int64)

    fast compute kurtosis series at each timewindow

# Input:
    - `v::RealArray`    : SeisData from SeisIO
    - `N::Int64`  : time window length to calculate kurtosis

    kurtosis evaluation following Baillard et al.(2013)
"""
function fast_kurtosis_series(v::Array, TN::Int64, SparseWin::Int64)

    kurt = zeros(length(v))
    n = length(v)
    kurt_grid = 1:n

    if n < TN error("Kurtosis time window is larger than data length. Decrease time window.") end
    if SparseWin > TN error("Sparse window is larger than Kurtosis time window. Decrease Kurtosis Sparse Window.") end

    #1. compute kurtosis with sparse grid
    kurt_sparse_grid = collect(TN:SparseWin:n)
    if kurt_sparse_grid[end] != n
        # fill up the edge of time series
        push!(kurt_sparse_grid, n)
    end

    kurt_sparse = []

    cm2 = 0.0  # empirical 2nd centered moment (variance)
    cm4 = 0.0  # empirical 4th centered moment

    t0 = @elapsed @simd for k = kurt_sparse_grid
        Trace = @views v[k-TN+1:k]
        m0 = mean(Trace)

        cm2 = Statistics.varm(Trace, m0, corrected=false)
        cm4 = fourthmoment(Trace, m0, corrected=false) #sum(xi - m)^4 / N
        #push!(kurt_sparse, (cm4 / (cm2 * cm2)) - 3.0)
        kurt[k] = (cm4 / (cm2 * cm2)) - 3.0
    end

    #2. interpolate kurtosis
    #t1 = @elapsed spl = Spline1D(kurt_sparse_grid, kurt_sparse; k=1, bc="nearest") #cubic spline
    #t2 = @elapsed kurt = evaluate(spl, kurt_grid)

    #println([t0, t1, t2])

    # !DEPRECATED!: 1. compute mean value at each time window by numerical sequence
    # 2. use mapreduce to sum values

    # # first term
    # Trace = @views v[1:TN]
    # z2 = zeros(TN)
    # m0 = mean(Trace)
    #
    # cm2 = Statistics.varm(Trace, m0, corrected=false)
    # cm4 = fourthmoment(Trace, m0, corrected=false)
    #
    # # fill first part with kurtosis at TN
    # kurt[1:TN] .= (cm4 / (cm2 * cm2)) - 3.0
    #
    # @simd for k = TN:n-1
    #
    #     diff1 = @inbounds @views (v[k-TN+1] - v[k+1])/TN
    #     m1 = m0 - diff1
    #     Trace = @views v[k-TN+2:k+1]
    #     cm2 = Statistics.varm(Trace, m1, corrected=false)
    #     cm4 = fourthmoment(Trace, m1, corrected=false) #sum(xi - m)^4 / N
    #     kurt[k+1] = (cm4 / (cm2 * cm2)) - 3.0
    #     m0 = m1
    # end

    # first term
    # Trace = @views v[1:TN]
    # z2 = zeros(TN)
    # m0 = mean(Trace)
    #
    # cm2 = Statistics.varm(Trace, m0, corrected=false)
    # cm4 = fourthmoment(Trace, m0, corrected=false)
    #
    # # fill first part with kurtosis at TN
    # kurt[1:TN] .= (cm4 / (cm2 * cm2)) - 3.0
    #
    # @simd for k = TN:n-1
    #
    #     diff1 = @inbounds @views (v[k-TN+1] - v[k+1])/TN
    #     m1 = m0 - diff1
    #     Trace = @views v[k-TN+2:k+1]
    #     cm2 = Statistics.varm(Trace, m1, corrected=false)
    #     cm4 = fourthmoment(Trace, m1, corrected=false) #sum(xi - m)^4 / N
    #     kurt[k+1] = (cm4 / (cm2 * cm2)) - 3.0
    #     m0 = m1
    # end

    return kurt

end


#---following functions are modified from Statistics.jl---#

centralizedabs4fun(m) = x -> abs2.(abs2.(x - m))
centralize_sumabs4(A::AbstractArray, m) =
    mapreduce(centralizedabs4fun(m), +, A)
centralize_sumabs4(A::AbstractArray, m, ifirst::Int, ilast::Int) =
    Base.mapreduce_impl(centralizedabs4fun(m), +, A, ifirst, ilast)

function centralize_sumabs4!(R::AbstractArray{S}, A::AbstractArray, means::AbstractArray) where S
    # following the implementation of _mapreducedim! at base/reducedim.jl
    lsiz = Base.check_reducedims(R,A)
    isempty(R) || fill!(R, zero(S))
    isempty(A) && return R

    if Base.has_fast_linear_indexing(A) && lsiz > 16 && !has_offset_axes(R, means)
        nslices = div(length(A), lsiz)
        ibase = first(LinearIndices(A))-1
        for i = 1:nslices
            @inbounds R[i] = centralize_sumabs4(A, means[i], ibase+1, ibase+lsiz)
            ibase += lsiz
        end
        return R
    end
    indsAt, indsRt = Base.safe_tail(axes(A)), Base.safe_tail(axes(R)) # handle d=1 manually
    keep, Idefault = Broadcast.shapeindexer(indsRt)
    if Base.reducedim1(R, A)
        i1 = first(Base.axes1(R))
        @inbounds for IA in CartesianIndices(indsAt)
            IR = Broadcast.newindex(IA, keep, Idefault)
            r = R[i1,IR]
            m = means[i1,IR]
            @simd for i in axes(A, 1)
                r += abs2(abs2(A[i,IA] - m))
            end
            R[i1,IR] = r
        end
    else
        @inbounds for IA in CartesianIndices(indsAt)
            IR = Broadcast.newindex(IA, keep, Idefault)
            @simd for i in axes(A, 1)
                R[i,IR] += abs2(abs2(A[i,IA] - means[i,IR]))
            end
        end
    end
    return R
end

function fourthmoment!(R::AbstractArray{S}, A::AbstractArray, m::AbstractArray; corrected::Bool=true) where S
    if isempty(A)
        fill!(R, convert(S, NaN))
    else
        rn = div(length(A), length(R)) - Int(corrected)
        centralize_sumabs4!(R, A, m)
        R .= R .* (1 // rn)
    end
    return R
end

"""
    fourthmoment(v, m; dims, corrected::Bool=true)
Compute the fourthmoment of a collection `v` with known mean(s) `m`,
optionally over the given dimensions. `m` may contain means for each dimension of
`v`. If `corrected` is `true`, then the sum is scaled with `n-1`,
whereas the sum is scaled with `n` if `corrected` is `false` where `n = length(v)`.
!!! note
    If array contains `NaN` or [`missing`](@ref) values, the result is also
    `NaN` or `missing` (`missing` takes precedence if array contains both).
    Use the [`skipmissing`](@ref) function to omit `missing` entries and compute the
    variance of non-missing values.
"""
fourthmoment(A::AbstractArray, m::AbstractArray; corrected::Bool=true, dims=:) = _fourthmoment(A, m, corrected, dims)

_fourthmoment(A::AbstractArray{T}, m, corrected::Bool, region) where {T} =
    fourthmoment!(Base.reducedim_init(t -> abs2(t)/2, +, A, region), A, m; corrected=corrected)

fourthmoment(A::AbstractArray, m; corrected::Bool=true) = _fourthmoment(A, m, corrected, :)

function _fourthmoment(A::AbstractArray{T}, m, corrected::Bool, ::Colon) where T
    n = length(A)
    n == 0 && return oftype((abs2(zero(T)) + abs2(zero(T)))/2, NaN)
    return centralize_sumabs4(A, m) / (n - Int(corrected))
end

end
