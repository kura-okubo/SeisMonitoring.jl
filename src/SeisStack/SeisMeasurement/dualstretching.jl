using Interpolations, Statistics, Distances, StatsBase, ColorSchemes, Plots
using SeisMonitoring: smooth_withfiltfilt
"""
    dualstretching(ref::AbstractArray, cur::AbstractArray, t::AbstractArray, fc::Float64,
                        window::AbstractArray;
                        dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500,
                        dQcinvmin::Float64=-0.01, dQcinvmax::Float64=+0.01, ntrial_q::Int=101,
                        dAAmin::Float64=-0.02, dAAmax::Float64=+0.02, ntrial_A::Int=101,
                        dist_method::String="euclidean")

Dual stretching method to estimate velocity and attenuation change.
Author: Kurama Okubo
"""
function dualstretching(ref::AbstractArray, cur::AbstractArray, t::AbstractArray, fc::Float64,
                        window::AbstractArray;
                        dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500,
                        dQcinvmin::Float64=-0.01, dQcinvmax::Float64=+0.01, ntrial_q::Int=101,
                        dAAmin::Float64=-0.02, dAAmax::Float64=+0.02, ntrial_A::Int=101,
                        dist_method::String="euclidean", figdir::String="", figname::String="", fillbox::AbstractArray=[])

    #---Process flow---#
    #1. compute dvv by stretching
    #2. compute envelope of ref and cur
    #3. compute dQinv and dss
    #------------------#

    #==NOTE:===#
    #Due to memory limitation, s_alltrace is available when !isempty(figdir) ==true (debug_plot=true)
    #===#

    MeasurementDict = Dict()

    #1. compute dvv by stretching
    dvv,cc_dvv,cdp,ϵ,allC = dualstretching_dvv(ref,cur,t,window,
                                    dvmin=dvmin, dvmax=dvmax,ntrial_v=ntrial_v)
    MeasurementDict["dvv"] = dvv
    MeasurementDict["cc_dvv"] = cc_dvv
    MeasurementDict["cdp_dvv"] = cdp
    MeasurementDict["ϵ_dvv"] = ϵ
    MeasurementDict["allC_dvv"] = allC

    #2. compute envelope of ref and cur
    Aref = abs.(hilbert(ref));
    Acur = abs.(hilbert(cur));

    #3. compute dss and dQinv
    dAA, dQcinv, dist_dualstretch, ϵA, ϵdQ, alldist_dualstretch, s_alltrace = dualstretching_dQinv(dvv, Aref, Acur, t, fc, window,
                dQcinvmin=dQcinvmin,
                dQcinvmax=dQcinvmax,
                dAAmin=dAAmin,
                dAAmax=dAAmax,
                ntrial_q=ntrial_q,
                ntrial_A=ntrial_A,
                dist_method = dist_method,
                figdir=figdir,
                figname=figname,
                fillbox=fillbox)

    MeasurementDict["dAA"] = dAA
    MeasurementDict["dQcinv"] = dQcinv
    MeasurementDict["dist_dualstretch"] = dist_dualstretch
    MeasurementDict["ϵA"] = ϵA
    MeasurementDict["ϵdQ"] = ϵdQ
    MeasurementDict["alldist_dualstretch"] = alldist_dualstretch
    # MeasurementDict["s_alltrace"] = s_alltrace # please avoid to save 3D array for the moment if the potential issue occurs in cluster.
    return MeasurementDict
end


function dualstretching_dvv(ref::AbstractArray,cur::AbstractArray,t::AbstractArray,
                    window::AbstractArray;
                    dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500)

    # NOTE: Instead of refining allC as written in https://github.com/lviens/2018_JGR,
    # test more stretching factor of current waveform using Cubic spline interpolation.

    ϵ = range(dvmin, stop=dvmax, length=ntrial_v)
    L = 1. .+ ϵ
    tau = t * L'
    allC = zeros(ntrial_v)

    # set of stretched/compressed current waveforms
    waveform_ref = ref[window]
    for ii = 1:ntrial_v
        # Note: stretching current waveform such that
        # u0(αt) = u1(t) (representation synchronized with this notebook)
        # compute correlation regime 1 -> regime 0
        #s = LinearInterpolation(tau[:,ii],cur,extrapolation_bc=Flat())(t)
        itp = Interpolations.interpolate(cur, BSpline(Cubic(Line(OnGrid()))))
        etpf = Interpolations.extrapolate(itp, Flat())
        tau_scale = range(tau[1,ii], step=abs(tau[2,ii]-tau[1,ii]), length=length(t))
        sitp = Interpolations.scale(etpf, tau_scale)
        s = sitp(t)
        waveform_cur = s[window]
        allC[ii] = Statistics.cor(waveform_ref,waveform_cur)
    end

    cdp = cor(cur[window],ref[window])
    # find the maximum correlation coefficient
    #dvv = 100. * ϵ[argmax(allC)]
    dvv = ϵ[argmax(allC)]
    cc = maximum(allC)

    return dvv,cc,cdp,Array(ϵ),allC
end

function dualstretching_dQinv(dvv::Float64, Aref::AbstractArray, Acur::AbstractArray, t::AbstractArray,fc::Float64,
            window::AbstractArray;
            dQcinvmin::AbstractFloat=-0.02,
            dQcinvmax::AbstractFloat=+0.02,
            dAAmin::AbstractFloat=-0.05,
            dAAmax::AbstractFloat=0.05,
            ntrial_q::Int=201,
            ntrial_A::Int=101,
            dist_method::String="euclidean",
            figdir::String="",
            figname::String="",
            fillbox::AbstractArray=[])

    stretch_debugplot = !isempty(figdir)
    ismahalanobis = (lowercase(dist_method) == "mahalanobis")

    ϵdQ = range(dQcinvmin, stop=dQcinvmax, length=ntrial_q)
    ϵA = range(dAAmin, stop=dAAmax, length=ntrial_A)
    α = (1.0 + dvv)
    β = (1.0 .+ ϵA)

    #---NOTE: 2020/06/11 Applying boxcar smoothing on Aref and Acur---#
    fs = 1.0/(t[2]-t[1])
    window_sec = 5.0 #[s]
    window_len = trunc(Int, fs*window_sec)
    Aref = smooth_withfiltfilt(Aref, window_len=window_len, window=:rect)
    Acur = smooth_withfiltfilt(Acur, window_len=window_len, window=:rect)
    #-------------------------------------------------------------#

    γ(β0, dQcinv0, fc, α0, t) = β0 * exp(-dQcinv0 * pi * fc * α0 * t)
    alldist = zeros(ntrial_A, ntrial_q)
    tau = α .* t # stretched time vector
    waveform_ref = Aref[window]

    # set of stretched/compressed current waveforms
    if stretch_debugplot
        s_alltrace = Array{Float64, 3}(undef, length(t), ntrial_A, ntrial_q)
    else
        s_alltrace = Array{Float64, 3}(undef, 0, 0, 0)
    end

    ismahalanobis && (Q_traces = Array{Float64, 2}(undef, length(window), ntrial_A*ntrial_q))

    itp = Interpolations.interpolate(Acur, BSpline(Cubic(Line(OnGrid()))))
    etpf = Interpolations.extrapolate(itp, Flat())
    tau_scale = range(tau[1], step=tau[2]-tau[1], length=length(t))
    sitp = Interpolations.scale(etpf, tau_scale)
    sitp_tr = sitp.(t)

    icount = 1

    t_trial = @elapsed for ii = 1:ntrial_A
        for jj = 1:ntrial_q
            # Note: stretching current envelope such that
            # γA0(αt, fc) = A1(t, fc)
            # then compute distance regime 1 -> regime 0

            #===NOTE: move the itp process outside of loop===#
            #=== as coda Qc stretching changes only amplitude===#
            # itp = Interpolations.interpolate(Acur, BSpline(Cubic(Line(OnGrid()))))
            # etpf = Interpolations.extrapolate(itp, Flat())
            # tau_scale = range(tau[1], step=tau[2]-tau[1], length=length(t))
            # sitp = Interpolations.scale(etpf, tau_scale)
            #================================================#

            s = zeros(length(t))
            for it = 1:length(t)
                s[it] = (1/γ(β[ii], ϵdQ[jj], fc, α, t[it])) * sitp_tr[it]
            end

            stretch_debugplot && (s_alltrace[:, ii, jj] = s)

            waveform_cur = s[window]

            if lowercase(dist_method) == "euclidean"
                alldist[ii, jj] = Distances.euclidean(waveform_ref, waveform_cur)
            elseif ismahalanobis
                Q_traces[:, icount] = waveform_cur
                icount += 1
            else
                error("dist method $(dist_method) is not available.")
            end
        end
    end

    if ismahalanobis
        #compute covariance matrix
        Q = StatsBase.cov(transpose(Q_traces))
        # println(size(Q))
        icount = 1
        for ii = 1:ntrial_A
            for jj = 1:ntrial_q
                alldist[ii, jj] = Distances.mahalanobis(waveform_ref, Q_traces[:,icount], Q)
                icount += 1
            end
        end
    end

    argmin(alldist)
    #cdp = cor(cur[window],ref[window])
    # find the maximum correlation coefficient
    dAA = ϵA[argmin(alldist)[1]]
    dQcinv = ϵdQ[argmin(alldist)[2]]
    dist_dualstretch = minimum(alldist)

    t_fig = @elapsed if !isempty(figdir)
        isempty(figname) || isempty(fillbox) && error("figname or fillbox is empty.")
        # plot stretching traces
        p_stretch = plot(bg=:white, dpi=100)

        # plot coda_window box
        yshift_box = [0,0]
        yrange = maximum(s_alltrace)
        plot!(fillbox[1:2], yshift_box,
            fillrange=[yshift_box.+yrange], fillalpha=0.1, c=:orange,
            label="", linealpha=0.0)
        plot!(fillbox[3:4], yshift_box,
            fillrange=[yshift_box.+yrange], fillalpha=0.1, c=:orange,
            label="", linealpha=0.0)

        plotspan = max(ceil(Int, size(s_alltrace, 2)/10), ceil(Int, size(s_alltrace, 3)/10))
        for ii = 1:plotspan:size(s_alltrace, 2) # loop in dAA
            for jj=1:plotspan:size(s_alltrace, 3) # loop in dQinv
                p_stretch = plot!(t, s_alltrace[:, ii, jj],
                color = ColorSchemes.rainbow[ii/size(s_alltrace, 2)],
                linealpha=jj/size(s_alltrace, 3), label="")
            end
        end

        # superplot reference
        p_stretch = plot!(t, Aref, color = :black, label="reference", legend=:topright)

        # superplot reference
        xlabel!(p_stretch, "Time lag [s]")
        ylabel!(p_stretch, "coherence")
        title!(p_stretch, "$(figname)")
        savefig(p_stretch, joinpath(figdir, "stretching_$(figname).png"))

        # plot heatmap of distance
        p_heat = heatmap(ϵdQ , ϵA, log10.(alldist), dpi=100)
        p_heat = Plots.xlabel!("Trial dQc^{-1}")
        p_heat = Plots.ylabel!("Trial dA/A")
        p_heat = Plots.title!("$(figname)")
        savefig(p_heat, joinpath(figdir, "heatmap_$(figname).png"))
    end

    # println("debug: t_trial = $(t_trial)[s], t_fig = $(t_fig)[s]")

    return dAA, dQcinv, dist_dualstretch, ϵA, ϵdQ, alldist, s_alltrace
end
#
#
# heatmap(dQcinv, ϵs, log10.(alldist))
# Plots.xlabel!("Estimated dQc^{-1}")
# Plots.ylabel!("Estimated ds/s")
# Plots.title!("Distance with examined dQinv and dss")
