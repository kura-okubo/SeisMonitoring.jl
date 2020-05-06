"""
    remove_eqfilt!(data::SeisChannel)

remove earthquake by kurtosis and STA/LTA threshold

# Input:
    - `data::SeisChannel`    : SeisData from SeisIO

"""
function remove_eqfilt!(data::SeisChannel, InputDict::OrderedDict)

    #
    # data_origin::SeisChannel, plot_kurtosis_α::Float64, max_taper_dur::Float64,
    # plot_boxheight::Float64, plot_span::Int64, plot_fmt::String, fodir::String, tstamp::String, tvec::Array{Float64,1}, IsSaveFig::Bool)

    noisesignal = data.misc["noisesignal"][:]
    nx = length(data.x)

    while i <= nx
        if !noisesignal[i]
            push!(t1, tvec[i])

            t1id = i

            #find next id
            tt1 += @elapsed nexttrueid = findfirst(x -> x == true, noisesignal[t1id:end])

            if isnothing(nexttrueid)
                # all data is removed
                t2id = nx
                iinc = nx
            else
                t2id = (t1id - 1) + nexttrueid - 1 # from first false index to end false index
                iinc = t2id - t1id + 1 # next true = end of false (t2id)+ 1
            end

            #=== Adaptive tukey windowing algorhitm below is written by Jared Bryan===#
            max_wintaper_length = Int(data.fs*InputDict["fixed_tukey_margin"])

            # compute α given maximum taper length
            eq_length = t2id-t1id+1
            taper_length = 2*max_wintaper_length
            tukey_length = eq_length + taper_length
            invert_tukey_α = taper_length/tukey_length

            # define inverted tukey window
            invtukeywin = -DSP.tukey(Int(tukey_length), invert_tukey_α) .+ 1.0

            # slice tukey window if it exceeds array bounds
            if t1id < max_wintaper_length + 1
                left_overflow = (max_wintaper_length-t1id)+1
                invtukeywin = @views invtukeywin[left_overflow+1:end]
                # leftmost t
                left = 1
            else
                # full taper length
                left = t1id - max_wintaper_length
            end

            if (nx - t2id + 1) < max_wintaper_length + 1
                right_overflow = (max_wintaper_length-(nx - t2id + 1))+1
                invtukeywin = @views invtukeywin[1:end-right_overflow]
                # rightmost t
                right = nx
            else
                # full taper length
                right = t2id + max_wintaper_length
            end
            #=================================================================#

            # apply tukey window
            data.x[left:right] .*= invtukeywin

        else
            iinc = 1
        end

        i += iinc

    end

end
