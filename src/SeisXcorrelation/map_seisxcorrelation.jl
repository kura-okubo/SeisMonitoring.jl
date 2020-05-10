
"""
map_seisxcorrelation(station_pair::String, InputDict::OrderedDict)

Compute cross-correlation save data in jld2 file with CorrData format.

# Arguments
- `data`::Dict    : !Reprecated! Dictionary containing input data in the form: tstamp/stn => SeisChannel
- `tstamp::String,`    : Time stamp to read from JLD2 file.
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters

# Output
- `tserrorList`::Array{String,1}    : Array containing the name of failed cross-correlation timestamp/stationpairs
- `basefoname.tstamp.jld2`    : contains CorrData structure with a hierarchical structure (CC function, metadata)

#
"""
function map_seisxcorrelation(key_station_pair::String, StationPairDict::OrderedDict, InputDict::OrderedDict)

    # open FileIO of RawData
    fi = jldopen(InputDict["cc_absolute_RawData_path"], "r")

    # open output file
    fopath = joinpath(InputDict["fodir"], key_station_pair*".jld2")
    #remove existing cc file
    ispath(fopath) && rm(fopath)
    fo = jldopen(fopath, "w")

    starts  = InputDict["starts"]
    ends    = InputDict["ends"]

    for channel_pair = StationPairDict[key_station_pair]
        sta1, sta2 = split(channel_pair, "-")
        # perform cross-correlation within the cc_time_unit window by window
        for tid in 1:length(starts)
            starttime, endtime = u2d.([starts[tid], ends[tid]])

            #1. assemble seisdata
            t_assemble = @elapsed S1, S2 = assemble_seisdata.([sta1, sta2], [fi], starttime, endtime, data_contents_fraction=InputDict["data_contents_fraction"])

            #2. applying phase shift so that the data is aligned at the sampling frequency
            phase_shift!.([S1, S2])

            #3. convert to RawData
            R1, R2 = RawData.([S1, S2], InputDict["cc_len"], InputDict["cc_step"])

            #4. detrend, taper and band pass before computing fft and cross-correlation
            clean_up!.([R1, R2], InputDict["freqency_band"][1], InputDict["freqency_band"][end])

            #5. apply one-bit normalization if true
            InputDict["IsOnebit"] && onebit!.([R1, R2])

            #NOTE: Future work; We can make an option here to move RawData to GPU:
            #if ["GPU"]; R1, R2 :> GPU; end

            #6. compute fft
            t_fft = @elapsed FFT1, FFT2 = compute_fft.([R1, R2])

            #7. spectral normalization
            InputDict["cc_method"] == "coherence" && coherence!.([FFT1, FFT2], InputDict["smoothing_half_win"], InputDict["waterlevel"])

            # first station is used as source; e.g. "BP.CCRB..BP1-BP.CCRB..BP1" then BP.CCRB..BP1 is used.
            InputDict["cc_method"] == "deconvolution" && deconvolution!(FFT1, InputDict["smoothing_half_win"], InputDict["waterlevel"])

            #8. Compute cross-correlation
            t_xcorr = @elapsed xcorr = compute_cc(FFT1, FFT2, InputDict["maxlag"], corr_type=InputDict["cc_method"])

            #9. Apply frequency-decomposion of cross-correlation function

            #10. Save corr-data





    close(fi)

end




















    # IO parameters
    finame     = InputDict["finame"]
    basefoname = InputDict["basefoname"]
    time_unit  = InputDict["timeunit"]
    # FFT parameters
    freqmin    = InputDict["freqmin"]
    freqmax    = InputDict["freqmax"]
    fs         = InputDict["fs"]
    cc_len     = InputDict["cc_len"]
    cc_step    = InputDict["cc_step"]
    to_whiten  = InputDict["to_whiten"]
    time_norm  = InputDict["time_norm"]
    half_win  = InputDict["half_win"]
    maxdistance= InputDict["maxdistance"]
    # correlation parameters
    corrtype   = InputDict["corrtype"]   # "xcorr," "acorr," or "xchancorr" (or any combination thereof)
    corrmethod = InputDict["corrmethod"] # "cross-correlation", "deconv", or "coherence"
    maxtimelag = InputDict["maxtimelag"]
    # stacking parameters
    stack      = InputDict["allstack"]
    max_std    = InputDict["max_std"]


    freqband = InputDict["freqband"]
    if typeof(freqband) == Int
        Nfreqband = freqband
    else
        Nfreqband = length(freqband)-1
    end

    # dictionary to cache FFTs
    FFTDict = Dict{String, FFTData}()
    # list of stations that failed for this timestep
    tserrorList = []

    # load input file
    infile = jldopen(finame, "r")

    # split dataset names (keys of data) by "/" to get station list
    # assume form "$tstamp/$station"
    #dsk = collect(keys(data))
    #stlist = sort([string(split.(dsk[i], "/")[2]) for i=1:length(dsk)])
    stlist =
    try
        collect(keys(infile["$tstamp"]))
    catch
        # this time stamp has no data
        # create empty output file for this time stamp, fill relevant info
        println("$tstamp has no data. skip this day with empty cc file $basefoname.$tstamp.jld2.")
        outfile = jldopen("$basefoname.$tstamp.jld2", "a+")
        outfile["info/stationlist"] = String[]
        outfile["info/timeunit"] = time_unit
        close(infile)
        close(outfile)
        return nothing
    end

    println("$tstamp: Computing cross-correlations")
    #stniter = 0 # counter to prevent computing duplicate xcorrs with reversed order
    # iterate over station list

    # clean preexisting cc file
    if ispath("$basefoname.$tstamp.jld2")
        rm("$basefoname.$tstamp.jld2")
    end

    outfile = jldopen("$basefoname.$tstamp.jld2", true, true, true, IOStream)
    outfile["info/stationlist"] = stlist
    outfile["info/timeunit"] = time_unit

    for (stniter, stn1) in enumerate(stlist)

        # don't attempt FFT if this failed already
        if stn1 in tserrorList continue end

        # try to read FFT from cached FFTs
        if stn1 in keys(FFTDict)
            FFT1 = FFTDict[stn1]
        else
            # read station SeisChannels into SeisData before FFT
            S1 = try
                    SeisData(infile["$tstamp/$stn1"])
                 catch
                     println("$tstamp: $stn1 encountered an error on FFT1 read seisdata. Skipping.")
                     push!(tserrorList, "$stn1")
                     continue
                 end

            if size(S1.t[1],1) > 2
                println("$tstamp: $stn1 has gap; for the moment, using ungap.")
                ungap!(S1, m=false)
                if length(S1[1].x) != S1.t[end, 1]
                    println("$tstamp: $stn1 has gap; ungap data size error. Skipping.")
                    push!(tserrorList, "$stn1")
                    continue
                end
            end

            if length(S1)[1] > 1 @warn "SeisData contains multiple channels. Operating only on the first." end
            delete!(S1[1].misc, "kurtosis")
            delete!(S1[1].misc, "eqtimewindow")

            # do not attempt fft if data was not available
            if haskey(S1[1].misc, "dlerror")
                if S1[1].misc["dlerror"] == 1
                    push!(tserrorList, "$stn1")
                    #filter!(a->a≠stn1, stlist)
                    println("$tstamp: $stn1 encountered an error: dlerror==1. Skipping.")
                    continue
                end
            end

            # make sure the data is the proper length to avoid dimension mismatch
            #npts1 = Int(time_unit * S1[1].fs)

            #if (length(S1[1].x) > npts1) S1[1].x=S1[1].x[1:npts1]; S1[1].t[end,1]=npts1 end
#--------
            phase_shift!(S1)
            R1 = RawData(S1,Int(cc_len), Int(cc_step))
            if time_norm == "one-bit"
                #println("apply one-bit normalization.")
                onebit!(R1)
            end
            #3. create FFTData

            if isnothing(R1)
                println("$tstamp: $stn1 encountered an error on FFT1: R1 nothing. Skipping.")
                push!(tserrorList, "$stn1")
                continue
            end

            tt1temp = @elapsed FFT1 = compute_fft(R1)
            # tt1temp = @elapsed FFT1 = compute_fft(S1, freqmin, freqmax, fs, Int(cc_step), Int(cc_len),
            #                                 to_whiten=to_whiten, time_norm=time_norm, max_std=max_std)
#--------
            # smooth FFT1 only if coherence is selected. Deconvolution will use only FFT2.
            if corrmethod == "coherence"
               coherence!(FFT1, half_win)
            end

            #println("fft1: $tt1temp ")
            FFTDict[stn1] = FFT1
        end

        # iterate over station list again
        for stn2 in stlist[stniter:end]
            if stn2 in tserrorList continue end # don't attempt FFT if this failed already

            # see if this is an auto-, cross-, or xchan-correlation
            ct = get_corrtype([stn1, stn2])

            # autocorrelation
            if (ct=="acorr-achan") && ("acorr-achan" in corrtype)
               # set the stn2 FFT to the already computed FFT for stn1
               FFT2 = deepcopy(FFT1)

            # cross- or cross-channel correlation
            elseif ct in corrtype

                # try to read FFT from cached FFTs
                if stn2 in keys(FFTDict)
                    FFT2 = FFTDict[stn2]
                else
                    # read station SeisChannels into SeisData before FFT
                    S2 = try
                            SeisData(infile["$tstamp/$stn2"])
                        catch
                            println("$tstamp: $stn2 encountered an error on FFT2 read seisdata. Skipping.")
                            push!(tserrorList, "$stn2")
                            continue
                        end

                    if size(S2.t[1],1) > 2
                        println("$tstamp: $stn2 has gap; for the moment, using ungap.")
                        ungap!(S2, m=false)
                        if length(S2[1].x) != S2.t[end, 1]
                            println("$tstamp: $stn2 has gap; ungap data size error. Skipping.")
                            push!(tserrorList, "$stn2")
                            continue
                        end
                    end

                    if length(S2)[1] > 1 @warn "SeisData contains multiple channels. Operating only on the first." end
                    delete!(S2[1].misc, "kurtosis")
                    delete!(S2[1].misc, "eqtimewindow")

                    # do not attempt fft if data was not available
                    if haskey(S2[1].misc, "dlerror")
                        if S2[1].misc["dlerror"] == 1
                            push!(tserrorList, "$stn2")
                            #filter!(a->a≠stn1, stlist)
                            println("$tstamp: $stn2 encountered an error: dlerror==1. Skipping.")
                            continue
                        end
                    end

                    # make sure the data is the proper length to avoid dimension mismatch
                    #npts2 = Int(time_unit * S2[1].fs)

                    #if (length(S2[1].x) > npts2) S2[1].x=S2[1].x[1:npts2]; S2[1].t[end,1]=npts2 end
#------------------------------------------------------------------#
                    # Phase shift SeisChannel if starttime is not aligned with sampling period.
                    phase_shift!(S2)
                    R2 = RawData(S2,Int(cc_len), Int(cc_step))
                    if time_norm == "one-bit"
                        #println("apply one-bit normalization.")
                        onebit!(R2)
                    end
                    #3. create FFTData

                    if isnothing(R2)
                        println("$tstamp: $stn2 encountered an error on FFT2: R2 nothing. Skipping.")
                        push!(tserrorList, "$stn2")
                        continue
                    end

                    tt2temp = @elapsed FFT2 = compute_fft(R2)
                    # tt2temp = @elapsed FFT2 = compute_fft(S2, freqmin, freqmax, fs, Int(cc_step), Int(cc_len),
                    #                 to_whiten=to_whiten, time_norm=time_norm, max_std=max_std)
#------------------------------------------------------------------#
                    # smooth FFT2 if it hasn't been smoothed already
                    if corrmethod == "coherence"
                        coherence!(FFT2, half_win)
                    end
                    #print("fft2: $tt2temp ")
                    FFTDict[stn2] = FFT2
                end
            else
                #println("Skipping cross-correlation of $stn1 and $stn2.")
                continue
            end

            # do not attempt fft if distance of station pair is too large (dist > maxdistance)
            #print(dist(S1[1].loc, S2[1].loc)./1e3)
            #println(" [km]")

            if maxdistance==false || dist(FFT1.loc, FFT2.loc) <= maxdistance

                # compute correlation using SeisNoise.jl -- returns type CorrData
                #println(FFT1.t)
                #println(FFT2.t)
                #println(intersect(FFT1.t,FFT2.t))
                # TODO: find a way to store FFT2 for deconvolution without storing two sets of FFTs (smooth and unsmoothed)
                # deconvolution divides the cross spectrum by the squared power spectrum of the second signal
                # but each signal can be both source and receiver
                if corrmethod == "deconv"
                    deconvolution!(FFT2, half_win)
                end

                tt3temp = @elapsed xcorr = compute_cc(FFT1, FFT2, maxtimelag, corr_type="cross-correlation")
                #print("xcorr: $tt3temp ")

		           #DEBUG
	            # if stn1 == stn2 == "BP.VCAB..BP1" && occursin("57", tstamp)
            	#    println(xcorr)
        		#    @show xcorr.corr
        		#    @show xcorr.misc
        		#    @show varname
        		#end

        		# try
                # compute distance between stations
                if isnothing(xcorr)
                    # compute_cc returns nothing, so skip this
                    continue;
                end

                #xcorr.misc["dist"] = dist(FFT1.loc, FFT2.loc)
                # save location of each station

                if stn1 == stn2
                    xcorr.misc["location"] = Dict(stn1=>FFT1.loc)
                else
                    xcorr.misc["location"] = Dict(stn1=>FFT1.loc, stn2=>FFT2.loc)
                end

                # stack over DL_time_unit
                if stack==true stack!(xcorr, allstack=true) end

                #DEBUG:
                # if occursin("2004.54", tstamp) && ct=="acorr"
                #     println("XCORR:")
                #     println(xcorr)
		        #     #continue
                # end

                # #===
                # 2020.02.13 Apply wavelet transform and append wtcorr to xcorr.misc
                # This inclease the strage use for cc, but will optimize cpu speed
                # 2020.02.26 move to just before dumping due to memory saving
                # ===#
    			# append_wtcorr!(xcorr, freqband, figdir="", α0 = InputDict["α0"], αmax = InputDict["αmax"])

                # push!(varnamelist, varname)
                # push!(xcorrlist, xcorr)

                # @show xcorrlist[end].misc["location"]

                # xcorr_of_varname = xcorrlist[i]

                # println("before append")
                # @show xcorr_of_varname.name
                # @show xcorr_of_varname.misc["location"]
                #DEBUG: xchan
                # debugcomp = xcorr.comp
                # println("cc name:"*xcorr.name)
                # println("cc comp:"*debugcomp*" must be = stn1:"*stn1*" stn2:"*stn2)
                # @show xcorr.misc["location"]

                # varname = "$tstamp/$stn1.$stn2"

                append_wtcorr!(xcorr, freqband, figdir="", α0 = InputDict["α0"], αmax = InputDict["αmax"])
                outfile["$tstamp/$stn1.$stn2"] = xcorr

                #
                # println("after append")
                # @show xcorr_of_varname.name
                # @show xcorr_of_varname.misc["location"]


                # jldopen("$basefoname.$tstamp.jld2", "a+") do outfile
                #     outfile[varname] = xcorr
                # end
                # catch y
                #     println(y)
                #     println("$stn1 and $stn2 have no overlap at $tstamp.")
                #     push!(tserrorList, "$stn1")
                # end
            end
        end

        # release memory held by the FFT and time series for this station
        delete!(FFTDict, stn1)
    end

    # save xcorr outside the loop
    outfile["info/errors"] = tserrorList

    close(infile)
    close(outfile)

    return nothing

end
