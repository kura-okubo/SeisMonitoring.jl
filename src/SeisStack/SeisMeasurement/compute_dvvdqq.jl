using Interpolations, Statistics, Distances, StatsBase, ColorSchemes, JLD2, DataFrames, StatsModels, GLM
using SeisMonitoring: smooth_withfiltfilt

"""
    compute_dvvdqq(ref::AbstractArray, cur::AbstractArray, t::AbstractArray, fc::Float64,
                        window::AbstractArray;
                        dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500,
                        dQcinvmin::Float64=-0.01, dQcinvmax::Float64=+0.01, ntrial_q::Int=101,
                        dAAmin::Float64=-0.02, dAAmax::Float64=+0.02, ntrial_A::Int=101,
                        dist_method::String="euclidean")

stretching and coda Q fitting method to estimate velocity and attenuation change.
Author: Kurama Okubo
"""
function compute_dvvdqq(ref::AbstractArray, cur::AbstractArray, t::AbstractArray, fc::Float64,
                        window::AbstractArray;
                        dvmin::Float64=-0.1,dvmax::Float64=0.1,ntrial_v::Int=500,
                        geometrical_spreading_α::Float64=0.5,
                        coda_smooth_window::Float64=10.0,
                        eps::Float64 = 1e-20,
                        figdir::String="",
                        figname::String="",
                        fillbox::AbstractArray=[])

    #---Process flow---#
    #1. compute dvv with stretching
    #2. compute dqq with linear regression
    #------------------#

    MeasurementDict = Dict()

    #1. compute dvv by stretching
    dvv,cc_dvv,cdp,ϵ,allC = stretching_dvv(ref,cur,t,window,
                                    dvmin=dvmin, dvmax=dvmax,ntrial_v=ntrial_v)
    MeasurementDict["dvv"] = dvv
    MeasurementDict["cc_dvv"] = cc_dvv
    MeasurementDict["cdp_dvv"] = cdp
    MeasurementDict["ϵ_dvv"] = ϵ
    MeasurementDict["allC_dvv"] = allC

    #2. compute dqq
    fs = 1.0/(t[2]-t[1]) # assuming uniform time space
    DqqDict = compute_dqq(dvv, ref, cur, t, fs, fc, geometrical_spreading_α,window,
                coda_smooth_window=coda_smooth_window,
                eps=eps,
                figdir=figdir,
                figname=figname,
                fillbox=fillbox)

    for key in keys(DqqDict)
        MeasurementDict[key] = DqqDict[key]
    end

    return MeasurementDict
end


function stretching_dvv(ref::AbstractArray,cur::AbstractArray,t::AbstractArray,
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
        allC[ii] = cor(waveform_ref,waveform_cur)
    end

    cdp = cor(cur[window],ref[window])
    # find the maximum correlation coefficient
    #dvv = 100. * ϵ[argmax(allC)]
    dvv = ϵ[argmax(allC)]
    cc = maximum(allC)

    return dvv,cc,cdp,Array(ϵ),allC
end


"""
    compute_dqq(dvv::Float64, tr_ref::AbstractArray, tr_cur::AbstractArray, t::AbstractArray, fs::Float64,
                fc::Float64,
                geometrical_spreading_α::Float64,
                coda_window::AbstractArray;
                coda_smooth_window::Float64=10.0,
                eps::Float64 = 1e-20,
                figdir::String="",
                figname::String="",
                fillbox::AbstractArray=[],)

compute dq/q and ds/s by exponential curve fitting

"""
function compute_dqq(dvv::Float64, tr_ref::AbstractArray, tr_cur::AbstractArray, t::AbstractArray, fs::Float64,
            fc::Float64,
            geometrical_spreading_α::Float64,
            coda_window::AbstractArray;
            coda_smooth_window::Float64=10.0,
            eps::Float64 = 1e-20,
            figdir::String="",
            figname::String="",
            fillbox::AbstractArray=[],)

    #1. Strech current trace by dvv, which is previously obtained using dv/v measurement.
    tau = t .* (1.0 + dvv)
    itp = Interpolations.interpolate(tr_cur, BSpline(Cubic(Line(OnGrid()))))
    etpf = Interpolations.extrapolate(itp, Flat())
    tau_scale = range(tau[1], step=tau[2]-tau[1], length=length(t))
    sitp = Interpolations.scale(etpf, tau_scale)
    dvvstretched_tr_cur = sitp.(t)

    #2. Compute reference and current Qc
    QcDict_ref = compute_codaQ(tr_ref, t, fs, geometrical_spreading_α, fc, coda_window, coda_smooth_window=coda_smooth_window, eps=eps)
    QcDict_cur = compute_codaQ(dvvstretched_tr_cur, t, fs, geometrical_spreading_α, fc, coda_window, coda_smooth_window=coda_smooth_window, eps=eps)

    #3. Compute dq/q
    # Refernce Qcinv

    Qcinv_pos_ref = QcDict_ref["Qcinv_pos"]
    Qcinv_neg_ref = QcDict_ref["Qcinv_neg"]

    # Current Qcinv
    Qcinv_pos_cur = QcDict_cur["Qcinv_pos"]
    Qcinv_neg_cur = QcDict_cur["Qcinv_neg"]

    # Compute dQcinv/Qcinv
    dqq_pos = (Qcinv_pos_cur - Qcinv_pos_ref)/Qcinv_pos_ref
    dqq_neg = (Qcinv_neg_cur - Qcinv_neg_ref)/Qcinv_neg_ref

    if !isnan(dqq_pos) && !isnan(dqq_neg)
        dqq_avg = (dqq_pos + dqq_neg)/2.0
    elseif !isnan(dqq_pos) && isnan(dqq_neg)
        dqq_avg = dqq_pos
    elseif isnan(dqq_pos) && !isnan(dqq_neg)
        dqq_avg = dqq_neg
    else
        dqq_avg = NaN
    end

    #4. Compute ds/s
    # Refernce S
    if !isnan(Qcinv_pos_ref)
        coef_pos_ref = coeftable(QcDict_ref["model_pos"]).cols[1]
        S_pos_ref = 10^(coef_pos_ref[1])
    else
        S_pos_ref = NaN
    end

    if !isnan(Qcinv_neg_ref)
        coef_neg_ref = coeftable(QcDict_ref["model_neg"]).cols[1]
        S_neg_ref = 10^(coef_neg_ref[1])
    else
        S_neg_ref = NaN
    end

    if !isnan(Qcinv_pos_cur)
        coef_pos_cur = coeftable(QcDict_cur["model_pos"]).cols[1]
        S_pos_cur = 10^(coef_pos_cur[1])
    else
        S_pos_cur = NaN
    end

    if !isnan(Qcinv_neg_cur)
        coef_neg_cur = coeftable(QcDict_cur["model_neg"]).cols[1]
        S_neg_cur = 10^(coef_neg_cur[1])
    else
        S_neg_cur = NaN
    end

    dss_pos = (S_pos_cur - S_pos_ref)/S_pos_ref
    dss_neg = (S_neg_cur - S_neg_ref)/S_neg_ref

    if !isnan(dss_pos) && !isnan(dss_neg)
        dss_avg = (dss_pos + dss_neg)/2.0
    elseif !isnan(dss_pos) && isnan(dss_neg)
        dss_avg = dss_pos
    elseif isnan(dss_pos) && !isnan(dss_neg)
        dss_avg = dss_neg
    else
        dss_avg = NaN
    end

    #5. Store parameter into dictionary
    DqqDict = Dict(
            "dqq_pos" => dqq_pos,
            "dqq_neg" => dqq_neg,
            "dqq_avg" => dqq_avg,
            "dss_pos" => dss_pos,
            "dss_neg" => dss_neg,
            "dss_avg" => dss_avg,
            "Qcinv_pos_ref" => Qcinv_pos_ref,
            "Qcinv_neg_ref" => Qcinv_neg_ref,
            "Qcinv_pos_cur" => Qcinv_pos_cur,
            "Qcinv_neg_cur" => Qcinv_neg_cur
            )

    t_fig = @elapsed if !isempty(figdir)

        !ispath(figdir) && mkpath(figdir)

        p1 = plot(bg=:white, size=(800, 400), dpi=100, legend=:topright)

        if !isempty(coda_window)

            # p1: plot raw trace and envelope function after correction of geometrical spreading
            yshift_box = [0,0]
            plot!(fillbox[1:2], yshift_box,
                fillrange=[yshift_box.-20], fillalpha=0.1, c=:orange,
                label="", linealpha=0.0)
            plot!(fillbox[1:2], yshift_box,
                fillrange=[yshift_box.+0.9], fillalpha=0.1, c=:orange,
                label="", linealpha=0.0)

            if length(fillbox) == 4
                plot!(fillbox[3:4], yshift_box,
                    fillrange=[yshift_box.-20], fillalpha=0.1, c=:orange,
                    label="", linealpha=0.0)
                plot!(fillbox[3:4], yshift_box,
                    fillrange=[yshift_box.+0.9], fillalpha=0.1, c=:orange,
                    label="", linealpha=0.0)
            end
        end

        # plot current abs cc, envelope and smoothed envelope
        # compute absolute (not hilbert envelope) signal energy
#         tr_curTa_log10 = log10.(abs.(tr_cur) .* abs.(t) .^geometrical_spreading_α)
        tr_curTa_log10 = log10.(abs.(dvvstretched_tr_cur) .* abs.(t) .^geometrical_spreading_α)

        plot!(t, tr_curTa_log10, label="abs(cc_cur*t^a)", color=:gray)
        plot!(t, QcDict_cur["Atα_log10"], label="Acur*t^a", color=:red, linestyle = :dash)
        plot!(t, QcDict_cur["Atα_log10_smoothed"], label="Acur*t^a smoothed", color=:red, linewidth=1.5)

        # plot reference smoothed envelope
        plot!(t, QcDict_ref["Atα_log10_smoothed"], label="Aref*t^a smoothed", color=:blue, linewidth=1.5)

        xmax = maximum(abs.(t))
        ymax = 0.95*maximum(QcDict_cur["Atα_log10"])
        plot!(xlim=(-xmax, xmax), ylim = (ymax-4, ymax), margin=5Plots.mm)
        title!("$(figname)")
        xlabel!("Time lag[s]")
        ylabel!("log10(Energy)")

        # p2: plot scaled fitting linear curve
        p2 = plot(bg=:white, size=(800, 400), dpi=100)

        if !isempty(coda_window)

            yshift_box = [0,0]
            plot!(fillbox[1:2], yshift_box,
                fillrange=[yshift_box.-20], fillalpha=0.1, c=:orange,
                label="", linealpha=0.0)
            plot!(fillbox[1:2], yshift_box,
                fillrange=[yshift_box.+0.9], fillalpha=0.1, c=:orange,
                label="", linealpha=0.0)

            if length(fillbox) == 4
                plot!(fillbox[3:4], yshift_box,
                    fillrange=[yshift_box.-20], fillalpha=0.1, c=:orange,
                    label="", linealpha=0.0)
                plot!(fillbox[3:4], yshift_box,
                    fillrange=[yshift_box.+0.9], fillalpha=0.1, c=:orange,
                    label="", linealpha=0.0)
            end
        end


        # plot linear fitting curve
        if !isnan(QcDict_ref["Qcinv_pos"])
            coef_pos_ref = coeftable(QcDict_ref["model_pos"]).cols[1]
            fit_curve_pos_ref = coef_pos_ref[2].* QcDict_ref["t_pos"]
            plot!(QcDict_ref["t_pos"], fit_curve_pos_ref, label="reference", color=:black)
        end

        if !isnan(QcDict_ref["Qcinv_meg"])
            coef_neg_ref = coeftable(QcDict_ref["model_neg"]).cols[1]
            fit_curve_neg_ref = coef_neg_ref[2].* QcDict_ref["t_neg"]
            plot!(QcDict_ref["t_neg"], fit_curve_neg_ref, label="reference", color=:black)
        end

        if !isnan(QcDict_cur["Qcinv_pos"])
            coef_pos_cur = coeftable(QcDict_cur["model_pos"]).cols[1]
            fit_curve_pos_cur = coef_pos_cur[2].* QcDict_cur["t_pos"]
            plot!(QcDict_cur["t_pos"], fit_curve_pos_cur, label="current", color=:red)
        end

        if !isnan(QcDict_cur["Qcinv_meg"])
            coef_neg_cur = coeftable(QcDict_cur["model_neg"]).cols[1]
            fit_curve_neg_cur = coef_neg_cur[2].* QcDict_cur["t_neg"]
            plot!(QcDict_cur["t_neg"], fit_curve_neg_cur, label="current", color=:red)
        end

        plot!(xlim=(-xmax, xmax), ylim = (-1.0, 0.5))
        xlabel!("Time lag[s]")
        ylabel!("log10(Energy)")
        title!("dqq_avg, dss_avg = $(dqq_avg), $(dss_avg)")
        p_all = plot(p1, p2, layout = (2, 1), size=(800, 800), margin=8Plots.mm)

        savefig(p_all, joinpath(figdir, "computedqq_$(figname).png"))

        # save QcDict_ref and QcDict_cur into jld2file
        jldopen(joinpath(figdir, "computedqq_$(figname).jld2"), "w") do file
            file["QcDict_ref"] = QcDict_ref
            file["QcDict_cur"] = QcDict_cur
        end
    end

    return DqqDict
end


"""
    compute_codaQ(x::AbstractArray, t::AbstractArray, fs::Float64, geometrical_spreading_α::Float64,
        fc::Float64, coda_window::AbstractArray, coda_smooth_window::Float64=10.0;eps::Float64 = 1e-20)

compute coda Q by exponential curve fitting

# Argument
- `x::AbstractArray`: cross_correlation time series
- `t::AbstractArray`: time lag
- `fs::Float64`     : Sampling Frequency
- `geometrical_spreading_α::Float64`: given geometrical spreading coefficient t^{-α}
- `fc::Float64`: mean frequency of the frequency band associated with cross_correlation function
- `coda_window::AbstractArray`: index of coda window to be fit
- `coda_smooth_window::AbstractArray`: [s] window length of boxcar smoothing
- `eps::Float64`: small number to avoid log10(0)

# Return
- `QcDict::Dict`: value of coda Q inverse and the other statistics.

"""
function compute_codaQ(x::AbstractArray, t::AbstractArray, fs::Float64, geometrical_spreading_α::Float64,
        fc::Float64, coda_window::AbstractArray;coda_smooth_window::Float64=10.0, eps::Float64 = 1e-20)

    #1. compute envelope*t^a
    A = abs.(hilbert(x));
    Atα = A .* abs.(t) .^geometrical_spreading_α

    #2. compute log 10
    Atα[Atα .< eps] .= eps # to avoid log10(0)
    Atα_log10 = log10.(Atα)

    # debug:this will be imported from SeisMonitoring
    # function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    #     w = getfield(DSP.Windows, window)(window_len)
    #     A = DSP.filtfilt(w ./ sum(w), A)
    #     return A
    # end

    #3. Apply boxcar smoothing
    window_len = trunc(Int, fs*coda_smooth_window)
    Atα_log10_smoothed = smooth_withfiltfilt(Atα_log10, window_len=window_len, window=:rect)

    #4. Split cc into positve and negative part
    t_pos, t_neg, A_pos, A_neg = split_cc(Atα_log10_smoothed, t)

    # tcenter = findfirst(x -> x >= 0.0, t)
    # pos_ind = tcenter:length(t)
    # neg_ind = tcenter:-1:1
    # t_pos   = t[pos_ind]
    # t_neg   = -t[neg_ind]
    # A_pos   = Atα_log10_smoothed[pos_ind]
    # A_neg   = Atα_log10_smoothed[neg_ind]

    # extract coda time window
    coda_pos_ind = findall(x -> x in t[coda_window], t_pos)
    coda_neg_ind = findall(x -> x in t[coda_window], t_neg)

    # make weights for linear regression
    # to consider antisymmetric coda window, process positive and negative side separately

    # Positive side
    if !isempty(coda_pos_ind)
        wts_pos = zeros(Float64, length(t_pos))
        wts_pos[coda_pos_ind] .= 1.0
        # linear regression using GLM module
        data_pos = DataFrame(X=t_pos, Y=A_pos)
        model_pos = GLM.lm(@formula(Y ~ X), data_pos, wts=wts_pos)
        coef_pos = coeftable(model_pos).cols[1]
        fit_curve_pos = coef_pos[1] .+ coef_pos[2].* t_pos
        #compute Qc inverse
        Qcinv_pos = (-coef_pos[2])/(pi * fc * log10(exp(1)))
    else
        wts_pos = zeros(Float64, length(t_pos))
        model_pos = []
        Qcinv_pos = NaN
    end


    if !isempty(coda_neg_ind)

        wts_neg = zeros(Float64, length(t_neg))
        wts_neg[coda_neg_ind] .= 1.0
        # linear regression using GLM module
        data_neg = DataFrame(X=t_neg, Y=A_neg)
        model_neg = GLM.lm(@formula(Y ~ X), data_neg, wts=wts_neg)
        coef_neg = coeftable(model_neg).cols[1]
        fit_curve_neg = coef_neg[1] .+ coef_neg[2].* t_neg
        #compute Qc inverse
        Qcinv_neg = (-coef_neg[2])/(pi * fc * log10(exp(1)))

    else
        wts_neg = zeros(Float64, length(t_neg))
        model_neg = []
        Qcinv_neg = NaN
    end

    if !isempty(coda_pos_ind) && !isempty(coda_neg_ind)
        Qcinv_avg = (Qcinv_pos+Qcinv_neg)/2
    elseif !isempty(coda_pos_ind) && isempty(coda_neg_ind)
        Qcinv_avg = Qcinv_pos
    elseif isempty(coda_pos_ind) && !isempty(coda_neg_ind)
        Qcinv_avg = Qcinv_neg
    else
        Qcinv_avg = NaN
    end

    QcDict = Dict("Qcinv_pos" => Qcinv_pos,
                "Qcinv_neg" => Qcinv_neg,
                "Qcinv_avg" => Qcinv_avg,
                "model_pos" => model_pos,
                "model_neg" => model_neg,
                "t_pos"=>t_pos,
                "t_neg"=>t_neg,
                "A_pos"=>A_pos,
                "A_neg"=>A_neg,
                "wts_pos"=>wts_pos,
                "wts_neg"=>wts_neg,
                "Atα_log10"=>Atα_log10,
                "Atα_log10_smoothed"=>Atα_log10_smoothed
                )

    return QcDict
end
