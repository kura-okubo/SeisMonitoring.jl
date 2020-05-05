__precompile__()
module Remove_eq

export detect_eq_kurtosis, stalta, remove_eq, s_whiten, s_whiten!

using SeisIO, SeisNoise, FFTW, DSP, StatsBase, Plots, JLD2, SeisIO, Printf

"""
    detect_eq_kurtosis(data::SeisChannel,tw::Float64=60.0, threshold::Float64=3.0, overlap::Float64=30)

find earthquake by kurtosis threshold

# Input:
    - `data::SeisChannel`    : SeisData from SeisIO
    - `tw::Float64`  : time window to evaluate if earthquake is contained.
    - `threshold::Float64` : kurtosis threshold: if kurtosis > threshold, the time window contains earthquake
    - `overlap::Float64`           : overlap of time window to control the margin of earthquake removal. (large overlap assigns large margin before and after earthquake.)

    kurtosis evaluation following Baillard et al.(2013)
"""
function detect_eq_kurtosis(data::SeisChannel; tw::Float64=60.0, kurtosis_threshold::Float64=3.0, overlap::Float64=30)

    #convert window lengths from seconds to samples
    twsize = trunc(Int, tw * data.fs)
    overlapsize = trunc(Int, overlap * data.fs)

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

        if !isnothing(findfirst(x -> x > kurtosis_threshold, twtrace))
            #this time window includes earthquake
            for tt= i+1:i+twsize
                data.misc["eqtimewindow"][tt] = false
            end
        end

        #advance long window
        i = i + slide

    end

    return data

end


"""
    detect_eq_stalta(data::SeisChannel,tw::Float64=60.0, threshold::Float64=3.0, overlap::Float64=30)

find earthquake and tremors by STA/LTA

# Input:
    - `data::SeisChannel`  : SeisChannel from SeisIO
    - `longwin::Float64`   : long time window
    - `shortwin::Float64`  : short time window
    - `threshold::Float64` : STA/LTA threshold: if STA/LTA > threshold, the time window contains earthquake
    - `overlap::Float64`   : overlap of time window to control the margin of earthquake removal. (large overlap assigns large margin before and after earthquake.)
    - `stalta_absoluteclip::Float64`   : clip if maximum absolute value exceeds this number (for the purpose of removing incoherent noise)

    original code written by Seth Olinger. For our purpose, overlap is applied for short time window
"""
function detect_eq_stalta(data::SeisChannel,longWinLength::Float64, shortWinLength::Float64, threshold::Float64, overlap::Float64, stalta_absoluteclip::Float64,
                            InputDict::Dict, tstamp::String, st1::String; datagap_eps::Float64=1e-8)


    #convert window lengths from seconds to samples
    longWin = trunc(Int,longWinLength * data.fs)
    shortWin = trunc(Int,shortWinLength * data.fs)
    overlapWin = trunc(Int,overlap * data.fs)
    trace = @view data.x[:]

    #buffer for sta/lta value
    staltatrace = zeros(length(trace))

    #calculate how much to move beginning of window each iteration
    slide = longWin

    #save weight
    eqweight = deepcopy(data.misc["eqtimewindow"])
    #print("before manipulation"); println(count(eqweight))

    # manipulate weight to avoid underestimation of lta due to data gap
    # threshold zero signal by 0.001% of lta within the timeseries, which encompasses small computational error such as 1e-21
    ltaall = StatsBase.mean(abs.(trace))
    for i = 1:length(trace)
        if isless(abs(trace[i]), ltaall*datagap_eps)
            # this is regarded as zero signal (data gap)
            eqweight[i] = false
        end
    end

    #print("after manipulation"); println(count(eqweight))

    #define long window counter and empty trigger vector
    i = 1

    #loop through current channel by sliding
    while i < length(trace)

        #check if last window and change long window length if so
        if length(trace) - i < longWin
            longWin = length(trace)-i
        end

        #define chunk of data based on long window length and calculate long-term average
        longTrace = trace[i:i+longWin]

        lta = StatsBase.mean(abs.(longTrace), weights(eqweight[i:i+longWin]))

        #reset short window counter
        n = 0
        breakflag = false
        #loop through long window in short windows
        while n <= longWin - shortWin

            #define chunk of data based on short window length and calculate short-term average
            shortTrace = @view trace[i+n:i+n+shortWin]
            shortTraceWeight = @view eqweight[i+n:i+n+shortWin]
            #sta = mean(abs.(shortTrace))
            sta = StatsBase.mean(abs.(shortTrace), weights(shortTraceWeight))
            #sta = StatsBase.mean(abs.(shortTrace))
            stamax = maximum(abs.(shortTrace))
            #calculate sta/lta ration
            staLta = sta/lta

            #record detection time if sta/lta ratio exceeds threshold
            if staLta > threshold || stamax > stalta_absoluteclip
                #triggers[i+n] = 1
                #this time window includes earthquake
                for tt= i+n:i+n+shortWin
                    data.misc["eqtimewindow"][tt] = false
                end
            end

            if breakflag
                break;
            end

            #advance short window
            n = n + shortWin - overlapWin

            #adjust for the last edge of long window
            if n >= longWin - shortWin
                n = longWin - shortWin
                breakflag = true
            end

            staltatrace[i+n+shortWin] = staLta

        end

        #advance long window
        i = i + slide

    end

    if InputDict["dumptraces"]
        #dump sta/lta trace
        mkpath(InputDict["dumppath"])
        tstamp_fname = replace(tstamp, ":" => "-")
        fname_out = join([tstamp_fname, st1, "stalta","dat"], '.')
        fopath = InputDict["dumppath"]*"/"*fname_out
        open(fopath, "w") do io
            for i in 1:length(staltatrace)
                write(io, @sprintf("%12.8f\n", staltatrace[i]))
            end
        end
    end
    return data

end


"""
    remove_eq(data::SeisChannel)

remove earthquake by kurtosis and STA/LTA threshold

# Input:
    - `data::SeisData`    : SeisData from SeisIO

"""
function remove_eq(data::SeisChannel, data_origin::SeisChannel, plot_kurtosis_α::Float64, max_taper_dur::Float64,
    plot_boxheight::Float64, plot_span::Int64, plot_fmt::String, fodir::String, tstamp::String, tvec::Array{Float64,1}, IsSaveFig::Bool)

    eqidlist = data.misc["eqtimewindow"][:]
    nx = length(data.x)

    i = 1

    t1 = []
    t2 = []
    y1 = []
    y2 = []

    tt1 = 0
    tt2 = 0
    tt3 = 0
    tt4 = 0
    tt5 = 0

    while i <= nx
        if !eqidlist[i]
            push!(t1, tvec[i])

            t1id = i

            #find next id
            tt1 += @elapsed nexttrueid = findfirst(x -> x == true, eqidlist[t1id:end])

            if isnothing(nexttrueid)
                # all data is removed
                t2id = nx
                iinc = nx
            else
                t2id = (t1id - 1) + nexttrueid - 1 # first false id to end false id
                iinc = t2id - i + 1 # first true as end of false + 1
            end

            push!(t2, tvec[t2id])

            max_wintaper_duration = Int(data.fs*max_taper_dur)

            # compute α given maximum taper length
            eq_length = t2id-t1id+1
            taper_length = 2*max_wintaper_duration
            tukey_length = eq_length + taper_length
            invert_tukey_α = taper_length/tukey_length

            # define inverted tukey window
            tt2 += @elapsed invtukeywin = -tukey(Int(tukey_length), invert_tukey_α) .+ 1

            # slice tukey window if it exceeds array bounds
            tt3 += @elapsed if t1id < max_wintaper_duration + 1
                left_overflow = (max_wintaper_duration-t1id)+1
                invtukeywin = @views invtukeywin[left_overflow+1:end]
                # leftmost t
                left = 1
            else
                # full taper length
                left = t1id - max_wintaper_duration
            end

            tt4 += @elapsed if (nx - t2id + 1) < max_wintaper_duration + 1
                right_overflow = (max_wintaper_duration-(nx - t2id + 1))+1
                invtukeywin = @views invtukeywin[1:end-right_overflow]
                # rightmost t
                right = nx
            else
                # full taper length
                right = t2id + max_wintaper_duration
            end

            # apply tukey window
            tt5 += @elapsed data.x[left:right] .*= invtukeywin

            #boxsize
            push!(y1, -plot_boxheight)
            push!(y2, plot_boxheight)

        else
            iinc = 1
        end

        i += iinc

    end

    #println([tt1, tt2, tt3, tt4, tt5])

    if IsSaveFig

        figdir = joinpath(fodir, "fig_removalEQ")
        if !ispath(figdir); mkpath(figdir); end

        Plots.plot(bg=:white)
        normalized_amp = 0.5 * maximum(data_origin.x[1:plot_span:end])

        Plots.plot!(tvec[1:plot_span:end], data_origin.x[1:plot_span:end]./ normalized_amp,
                    label="raw data", color="black")

        Plots.plot!(tvec[1:plot_span:end], data.x[1:plot_span:end]./ normalized_amp,
                    label="after removal", color="blue")

        #to plot kurtosis computed points
        kurtid = findall(x -> !iszero(data.misc["kurtosis"][x]), 1:length(data.misc["kurtosis"]))
        if isempty(filter(!isnan, data.misc["kurtosis"][kurtid]))
            kurtnormalize = 1.0
        else
            kurtnormalize = maximum(filter(!isnan, data.misc["kurtosis"][kurtid]))
        end

        Plots.plot!(tvec[kurtid], data.misc["kurtosis"][kurtid] ./ kurtnormalize .* plot_kurtosis_α,
                    color="red",
                    label="kurtosis")

        p = Plots.plot!(xlabel = "[hours]",
                    ylabel = "[m/s]",
                    title =  @sprintf("%s %s", data.id, tstamp),
                    size = (1200, 800))

        figname = @sprintf("%s/%s_%s.%s", figdir, data.id, tstamp, plot_fmt)
        Plots.savefig(p, figname)


        # save jld2 fig file
        # trace1 = scatter(;x=tvec[1:plot_span:end], y=data_origin.x[1:plot_span:end]./ normalized_amp,
        #  mode="lines", name="raw data", line_color="black")
        #
        # trace2 = scatter(;x=tvec[1:plot_span:end], y=data.x[1:plot_span:end]./ normalized_amp,
        #  mode="lines", name="after remove", line_color="blue")
        #
        # #to plot kurtosis computed points
        # kurtid = findall(x -> !iszero(data.misc["kurtosis"][x]), 1:length(data.misc["kurtosis"]))
        #
        # if isempty(filter(!isnan, data.misc["kurtosis"][kurtid]))
        #     kurtnormalize = 1.0
        # else
        #     kurtnormalize = maximum(filter(!isnan, data.misc["kurtosis"][kurtid]))
        # end
        #
        # trace3 = PlotlyJS.scatter(;x=tvec[kurtid], y=data.misc["kurtosis"][kurtid] ./ kurtnormalize .* plot_kurtosis_α,
        #  line_color="red", mode="lines", name="kurtosis")
        #
        # if !isempty(t1)
        #  shapes = PlotlyJS.rect(t1, t2, y1, y2; fillcolor="#ff99cc", opacity=0.3, line_width=0)
        #  layout = PlotlyJS.Layout(shapes=shapes, width=1200, height=600,
        #      xaxis=attr(title="Time [hour]"),
        #      yaxis=attr(title="Normalized velocity"),
        #      font =attr(size=13),
        #      showlegend=true,
        #      title = @sprintf("%s %s", data.id, tstamp))
        #
        #      p = PlotlyJS.plot([trace1; trace2; trace3],layout)
        #
        # else
        #  layout = PlotlyJS.Layout(width=1200, height=600,
        #      xaxis=attr(title="Time [hour]"),
        #      yaxis=attr(title="Normalized velocity"),
        #      font =attr(size=12),
        #      showlegend=true,
        #      title = @sprintf("%s %s", data.id, tstamp))
        #
        #  p = PlotlyJS.plot([trace1; trace2; trace3],layout)
        #
        # end
        #
        # figdir = joinpath(fodir, "fig")
        # mkpath(figdir)
        # figname = @sprintf("%s/%s_%s.%s", figdir, data.id, tstamp, plot_fmt)
        # ORCA.savefig(p, figname)
        # #display(p)
        # #println("press return for next plot...")
        # #readline()
    end

    return data

end

"""
    s_whiten!(data::SeisChannel)

Apply spectral whitening to seischannel trace
# Input:
    - `data::SeisChannel`    : SeisChannel from SeisIO
    - `freqmin::Real`: Pass band low corner frequency.
    - `freqmax::Real`: Pass band high corner frequency.
    - `pad::Int`: Number of tapering points outside whitening band.
"""
function s_whiten!(data::SeisChannel, freqmin::Float64, freqmax::Float64;pad::Int=50)

    # compute fft of time series
    FFT = rfft(data.x,1)

    # to use SeisNoise.whiten!() prepare (N, 2) array
    FFTtemp = Array{Complex{Float32}}(undef, length(FFT), 2)
    FFTtemp[:, 1] = FFT

    SeisNoise.whiten!(FFTtemp,freqmin,freqmax,data.fs, data.t[end, 1], pad=pad)

    data.x = irfft(FFTtemp[:,1],data.t[end, 1],1)
    return nothing

end
s_whiten(data::SeisChannel, freqmin::Float64, freqmax::Float64;pad::Int=50) = (U = deepcopy(data);
    s_whiten!(U,freqmin,freqmax,pad=pad);
    return U)

end
