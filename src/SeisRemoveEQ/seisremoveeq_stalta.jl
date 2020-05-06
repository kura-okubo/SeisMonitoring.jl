"""
    detect_eq_stalta(data::SeisChannel, InputDict::OrderedDict)

find earthquake and tremors using STA/LTA
original code written by Seth Olinger. For our purpose, overlap is applied to short time window
"""
function detect_eq_stalta(data::SeisChannel, InputDict::OrderedDict, datagap_eps::Float64=1e-8)

    #
    # longWinLength::Float64, shortWinLength::Float64, threshold::Float64, overlap::Float64, stalta_absoluteclip::Float64,
    #                         InputDict::Dict, tstamp::String, st1::String; datagap_eps::Float64=1e-8)


    longWinLength       = InputDict["stalta_longwindow"]
    shortWinLength      = InputDict["stalta_shortwindow"]
    stalta_threshold    = InputDict["stalta_threshold"]
    stalta_overlap      = InputDict["stalta_overlap"]
    stalta_absoluteclip = InputDict["stalta_absoluteclip"]

    #convert window lengths from seconds to samples
    longWin      = trunc(Int,longWinLength * data.fs)
    shortWin     = trunc(Int,shortWinLength * data.fs)
    overlapWin   = trunc(Int,overlap * data.fs)
    trace        = @view data.x[:]

    #buffer for sta/lta value
    staltatrace = zeros(length(trace))

    #calculate how much to move beginning of window each iteration
    slide = longWin

    #save weight
    nonzero_eqweight = deepcopy(data.misc["noisesignal"])

    #NOTE: manipulate weight to avoid underestimation of lta due to data gap
    # threshold zero signal by datagap_eps::Float64=1e-8 of lta within the timeseries,
    # which encompasses small computational error such as 1e-21
    ltaall = StatsBase.mean(abs.(trace))
    for i = 1:length(trace)
        if isless(abs(trace[i]), ltaall*datagap_eps)
            # this is regarded as zero signal (data gap)
            nonzero_eqweight[i] = false
        end
    end

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

        lta = StatsBase.mean(abs.(longTrace), weights(nonzero_eqweight[i:i+longWin]))

        #reset short window counter
        n = 0
        breakflag = false
        #loop through long window in short windows
        while n <= longWin - shortWin

            #define chunk of data based on short window length and calculate short-term average
            shortTrace = @view trace[i+n:i+n+shortWin]
            shortTraceWeight = @view nonzero_eqweight[i+n:i+n+shortWin]
            #sta = mean(abs.(shortTrace))
            sta = StatsBase.mean(abs.(shortTrace), weights(shortTraceWeight)) # weights is applied to
            #sta = StatsBase.mean(abs.(shortTrace))
            stamax = maximum(abs.(shortTrace))
            #calculate sta/lta ration
            staLta = sta/lta

            #record detection time if sta/lta ratio exceeds threshold
            if staLta > threshold || stamax > stalta_absoluteclip
                #this time window includes earthquake
                for tt= i+n:i+n+shortWin
                    data.misc["noisesignal"][tt] = false
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

    #append stalta trace
    data["stalta_trace"] = staltatrace

end
