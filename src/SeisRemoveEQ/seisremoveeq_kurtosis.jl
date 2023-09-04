
"""
    get_kurtosis(data::SeisChannel,kurtsis_tw_sparse::Float64; timewinlength::Float64=60)

    compute kurtosis at each timewindow

# Input:
    - `data::SeisChannel`    : SeisData from SeisIO
    - `kurtosis_tw_sparse::Float64` : time length of span for kurtosis time window
    - `timewinlength::Float64`  : time window to calculate kurtosis
    kurtosis evaluation following Baillard et al.(2013)
"""
function get_kurtosis!(data::SeisChannel, InputDict::OrderedDict)

    timewinlength = InputDict["shorttime_window"]
    kurtosis_tw_sparse = InputDict["timewindow_overlap"]

    #convert window lengths from seconds to samples
    TimeWin = trunc(Int, timewinlength * data.fs)
    SparseWin = trunc(Int, kurtosis_tw_sparse * data.fs)
    data.misc["kurtosis"] = fast_kurtosis_series(data.x, TimeWin, SparseWin)

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

    if n < TN; error("Kurtosis time window is larger than data length. Decrease time window.") end
    if SparseWin > TN; error("Sparse window is larger than Kurtosis time window. Decrease Kurtosis Sparse Window.") end

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

    return kurt

end


#---following functions are modified from Statistics.jl (https://github.com/JuliaLang/Statistics.jl)---#

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

#----------------------------------------------------------------------------#

"""
    detect_eq_kurtosis!(data::SeisChannel,tw::Float64=60.0, threshold::Float64=3.0, overlap::Float64=30)

find earthquake by kurtosis threshold

# Input:
    - `data::SeisChannel`    : SeisData from SeisIO
    - `tw::Float64`  : time window to evaluate if earthquake is contained.
    - `threshold::Float64` : kurtosis threshold: if kurtosis > threshold, the time window contains earthquake
    - `overlap::Float64`           : overlap of time window to control the margin of earthquake removal. (large overlap assigns large margin before and after earthquake.)

    kurtosis evaluation following Baillard et al.(2013)
    Algorithm for the loop through current channel by slidingis is contributed by Seth Olinger.
"""
function detect_eq_kurtosis!(data::SeisChannel, InputDict::OrderedDict)

    kurtosis_removewindow   = InputDict["shorttime_window"]
    kurtosis_threshold      = InputDict["kurtosis_threshold"]
    kurtosis_overlap        = InputDict["timewindow_overlap"]

    #convert window lengths from seconds to samples
    twsize = trunc(Int, kurtosis_removewindow * data.fs)
    overlapsize = trunc(Int, kurtosis_overlap * data.fs)

    if twsize == overlapsize
        error("shortWinLength should be longer than stalta_overlap.")
    end

    #calculate how much to move beginning of window each iteration
    slide = twsize-overlapsize

    #kurtosis of timeseries
    ku1 = data.misc["kurtosis"][:]
    #reset long window counter and triggers for current channel
    i = 0

    #loop through current channel by sliding
    while i < length(ku1) - twsize

        #check if last window and change long window length if so
        if length(ku1) - i < twsize
            twsize = length(ku1)-i
        end

        #define chunk of data based on long window length and calculate long-term average
        twtrace = @views ku1[i+1:i+twsize]

        if any(x -> x > kurtosis_threshold, twtrace)
            #this time window includes earthquake
            for tt= i+1:i+twsize
                data.misc["noisesignal"][tt] = false
            end
        end

        #advance long window
        i = i + slide

    end

end
