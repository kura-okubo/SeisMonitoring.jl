"""
    get_noisedatafraction(A::AbstractArray; zerosignal_minpts=100, eps_α=1e-6)

Get fraction of noise contents in waveform.
Zero signal is detected if the signal with its absolute value is less than eps_α * mean(A)
continues more than zerosignal_minpts.

# Return
- `datacontentsfraction::Float64`: fraction of noise data containts with the trace.
"""
function get_noisedatafraction(A::AbstractArray; zerosignal_minpts=100, eps_α=1e-6)

    signal_threshold = eps_α * mean(abs.(A))
    nullsignals = findall(x -> abs(x) < signal_threshold, A)

    # find nullsignals which is continuous more than zerosignal_minpts.
    i = 1

    Nnullsignal = 0
    seriescount = 1
    @simd for i = 1:length(nullsignals)-1
        # println(seriescount)
        # if nullsignals[i] is continuous to next index
        if nullsignals[i] + 1 == nullsignals[i+1]
            # this is continuous
            seriescount += 1

            # manipulate the last series
            if i==length(nullsignals)-1
                seriescount >= zerosignal_minpts && (Nnullsignal += seriescount)
            end

        else
            # append number of nullsignals if it is more than zerosignal_minpts, then reset series count
            seriescount >= zerosignal_minpts && (Nnullsignal += seriescount)
            seriescount = 1
        end
    end

    return 1.0 - (float(Nnullsignal) / float(length(A)))

end

get_noisedatafraction(Ch::SeisChannel;zerosignal_minpts=100, eps_α=1e-6) =
(return get_noisedatafraction(Ch.x, zerosignal_minpts=zerosignal_minpts, eps_α=eps_α))
