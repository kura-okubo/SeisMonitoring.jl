#
# """
# defaultinputdict(InputDict::Dict)
#
# default if input parameter is missing.
# """
# function defaultinputdict!(InputDict::Dict)
#
# 	# default values
# 	def = Dict()
# 	def["finame"]				 	= "./input.jld2"
# 	def["IsKurtosisRemoval"] 		= true
# 	def["max_edgetaper_duration"] 	= 60 * 5
# 	def["kurtosis_tw_sparse"] 		= 60
# 	def["kurtosis_timewindow"] 		= 60*3
# 	def["kurtosis_threshold"] 		= 2.0
# 	def["IsSTALTARemoval"] 			= true
# 	def["stalta_longtimewindow"] 	= 60*60*2
# 	def["stalta_threshold"] 		= 1.2
# 	def["stalta_absoluteclip"] 		= 1e20
# 	def["max_wintaper_duration"] 	= 60 * 3
# 	def["removal_shorttimewindow"] 	= 60 * 3
# 	def["overlap"] 					= 60
# 	def["IsIsolateComponents"] 		= false
# 	def["priority_channles"] 		= Dict()
# 	def["IsSaveFig"] 				= false
# 	def["plot_kurtosis_α"] 			= 1.2
# 	def["plot_boxheight	"] 			= 1.5
# 	def["plot_span"] 				= 100
# 	def["outputformat"]				= "JLD2"
# 	def["IsStartendtime"] 			= false
# 	def["fodir"] 					= "./dataset"
# 	def["foname"] 					= "eq_removed.jld2"
#
# 	#For dump everytrace
# 	#This is always false. If you need to plot some figures for
# 	#details of kurtosis and STA/LTA, please turn it true; it slowdowns down computation.
# 	def["dumptraces"] 				= false
# 	def["dumppath"]					= InputDict["fodir"]*"/dumptraces"
#
# 	for key in keys(def)
# 		if !haskey(InputDict, key)
# 			InputDict["$key"] = def["$key"]
# 		end
# 	end
#
# end

# """
# initlogo()
#
# print initial logo
# """
# function initlogo()
#
#     print("
#     _____      _      ______
#    /  ___|    (_)     | ___ \\
#    \\ `--.  ___ _ ___  | |_/ /___ _ __ ___   _____   _____
#     `--. \\/ _ \\ / __| |    // _ \\ '_ ` _ \\ / _ \\ \\ / / _ \\
#    /\\__/ /  __/ \\__ \\ | |\\ \\  __/ | | | | | (_) \\ V /  __/
#    \\____/ \\___|_|___/ \\_| \\_\\___|_| |_| |_|\\___/ \\_/ \\___|
#     ______ ____
#    |  ____/ __ \\       |
#    | |__ | |  | |      |
#    |  __|| |  | |      |  v1.0 (Last update 07/07/2019)
#    | |___| |__| |      |  © Kurama Okubo
#    |______\\___\\_\\      |
#
# ")
#
#     println("Job start running at "*string(now())*"\n")
#
# end
#
# """
# loaddata(finame, path)
# """
# function loaddata(finame, path)
#     try
#         @suppress_err return FileIO.load(finame, path)
#     catch
#         return false;
#     end
# end



        # finame                 =InputDict["finame"]
        # IsKurtosisRemoval      =InputDict["IsKurtosisRemoval"]
        # max_edgetaper_duration =InputDict["max_edgetaper_duration"]
        # kurtosis_tw_sparse     =InputDict["kurtosis_tw_sparse"]
        # kurtosis_timewindow    =InputDict["kurtosis_timewindow"]
        # kurtosis_threshold     =InputDict["kurtosis_threshold"]
        # IsSTALTARemoval        =InputDict["IsSTALTARemoval"]
        # stalta_longtimewindow  =InputDict["stalta_longtimewindow"]
        # stalta_threshold       =InputDict["stalta_threshold"]
        # stalta_absoluteclip    =InputDict["stalta_absoluteclip"]
        # max_wintaper_duration  =InputDict["max_wintaper_duration"]
        # removal_shorttimewindow=InputDict["removal_shorttimewindow"]
        # overlap                =InputDict["overlap"]
        # Iswhiten               =InputDict["Iswhiten"]
        # freqmin_whiten         =InputDict["freqmin_whiten"]
        # freqmax_whiten         =InputDict["freqmax_whiten"]
        # plot_kurtosis_α        =InputDict["plot_kurtosis_α"]
        # plot_boxheight         =InputDict["plot_boxheight"]
        # plot_span              =InputDict["plot_span"]
        # plot_fmt               =InputDict["plot_fmt"]
        # fodir                  =InputDict["fodir"]
        # foname                 =InputDict["foname"]
        # fopath                 =InputDict["fopath"]
        # IsSaveFig              =InputDict["IsSaveFig"]
        #
        #
        # DLtimestamplist        =InputDict["DLtimestamplist"]
        # stationlist            =InputDict["stationlist"]
        # NumofTimestamp         =InputDict["NumofTimestamp"]
        #
        # tstamp = DLtimestamplist[dlid]
        #
        # if mod(dlid, round(0.1*NumofTimestamp)+1) == 0
        #     println(@sprintf("start process %s", tstamp))
        # end
        #
        # bt_getkurtosis = 0.0
        # bt_removeeq = 0.0
        #
        # for st = stationlist
        #     #S = t[joinpath(tstamp, st)]
        #     st1 = replace(st, "-"=>"")
        #
        #     # load raw data
        #     Stemp = loaddata(finame, joinpath(tstamp, st1))
        #
        #     # skip if it does not exist
        #     if Stemp == false
        #         continue;
        #     end
        #
        #     #convert if it's SeisChannel
        #     if Stemp isa SeisIO.SeisChannel
        #         Stemp = SeisData(Stemp)
        #     end
        #
        #     SremEQ = SeisData()
        #
        #     #loop by SeisChannel and save into temporal wile
        #     for Seisid = 1:Stemp.n
        #
        #         Schan = Stemp[Seisid]
        #
        #         if !haskey(Schan.misc, "dlerror")
        #             if iszero(Schan.x)
        #                 Schan.misc["dlerror"] = true
        #             else
        #                 Schan.misc["dlerror"] = false
        #             end
        #         end
        #
        #         #---Transient error on combination between download margin and resampling---#
        #         # S.t may contain index 0 in t[:, 1]. in that case , this causes the error on later part, so discard that data.
        #         if any(Schan.t[:,1] .== 0)
        #             println("data gap index found.")
        #             println(Schan)
        #             Schan.misc["dlerror"] = true
        #         end
        #         #---------------------------------------------------------------------------#
        #
        #         if Schan.misc["dlerror"] == 0
        #
        #             dt = 1/Schan.fs
        #             tvec = collect(0:length(Schan.x)-1) * dt ./ 60 ./ 60
        #
        #             #tapering to avoid instrumental edge artifacts
        #             #SeisIO.taper!(Schan,  t_max = max_edgetaper_duration, α=0.05) this is already done by download margin
        #
        #             S1 = deepcopy(Schan)
        #
        #             #set long window length to user input since last window of previous channel will have been adjusted
        #             S1.misc["eqtimewindow"] = fill(true, length(S1.x))
        #
        #             if IsKurtosisRemoval
        #                 # compute kurtosis and detect
        #                 bt_1 = @elapsed S1 = Get_kurtosis.get_kurtosis(S1, float(kurtosis_timewindow), float(kurtosis_tw_sparse))
        #
        #                 bt_2 = @elapsed S1 = Remove_eq.detect_eq_kurtosis(S1, tw=float(removal_shorttimewindow), kurtosis_threshold=float(kurtosis_threshold), overlap=float(overlap))
        #
        #                 btsta_1 = 0
        #
        #                 if IsSTALTARemoval
        #                     # detect earthquake and tremors by STA/LTA
        #                     btsta_1 = @elapsed S1 = Remove_eq.detect_eq_stalta(S1, float(stalta_longtimewindow), float(removal_shorttimewindow),
        #                                         float(stalta_threshold), float(overlap), float(stalta_absoluteclip),  InputDict, tstamp, st1)
        #                 end
        #
        #                 if Iswhiten
        #                     Remove_eq.s_whiten!(S1, freqmin_whiten, freqmax_whiten)# apply spectral whitening on this channel
        #                 end
        #
        #                 bt_3 = @elapsed S1 = Remove_eq.remove_eq(S1, Schan, float(plot_kurtosis_α), float(max_wintaper_duration),
        #                                 plot_boxheight, trunc(Int, plot_span), plot_fmt, fodir, tstamp, tvec, IsSaveFig)
        #
        #                 bt_getkurtosis += bt_1
        #                 bt_removeeq += bt_2 + bt_3 + btsta_1
        #
        #                 #println([bt_2, bt_3, btsta_1])
        #                 #if mod(dlid, round(0.1*NumofTimestamp)+1) == 0
        #                 #    println([bt_1, bt_2, btsta_1, bt_3])
        #                 #end
        #
        #             else
        #                 #only STA/LTA
        #                 if IsSTALTARemoval
        #
        #                     data.misc["kurtosis"] = zeros(length(S1.x))
        #
        #                     bt_2 = @elapsed S1 = Remove_eq.detect_eq_stalta(S1, float(stalta_longtimewindow), float(removal_shorttimewindow),
        #                                         float(stalta_threshold), float(overlap), InputDict, tstamp, st1)
        #
        #                     if Iswhiten
        #                         Remove_eq.s_whiten!(S1, freqmin_whiten, freqmax_whiten) # apply spectral whitening on this channel
        #                     end
        #
        #                     bt_3 = @elapsed S1 = Remove_eq.remove_eq(S1, S, float(plot_kurtosis_α), float(max_wintaper_duration),
        #                                     plot_boxheight, trunc(Int, plot_span), plot_fmt, fodir, tstamp, tvec, IsSaveFig)
        #
        #
        #                     bt_getkurtosis += 0.0
        #                     bt_removeeq += bt_2 + bt_3
        #
        #                 else
        #                     #no removal process.
        #                     @warn "Both 'IsKurtosisRemoval' and 'IsKurtosisRemoval' are false. No removal process is executed. Abort."
        #                     exit(0)
        #                 end
        #
        #             end
        #
        #             if InputDict["IsOutputRemovalFrac"]
        #                 #output removal fraction on this channel
        #                 eqidlist = S1.misc["eqtimewindow"][:]
        #                 numofremoval = sum(x->x==false, eqidlist, dims=1)
        #                 fractionofremoval = numofremoval[1] / length(eqidlist)
        #
        #                 y, jd = parse.(Int64, split(InputDict["DLtimestamplist"][dlid], ".")[1:2])
        #                 tstamp_fname = replace(tstamp, ":" => "-")
        #                 fname_out = join([tstamp_fname, st1, "removalfraction","dat"], '.')
        #                 fopath = InputDict["removal_fractionpath"]*"/"*fname_out
        #                 open(fopath, "w") do io
        #                    write(io, @sprintf("%f\n", fractionofremoval))
        #                 end
        #             end
        #
        #
        #             if InputDict["dumptraces"]
        #                 #dump raw trace
        #                 tstamp_fname = replace(tstamp, ":" => "-")
        #                 fname_out = join([tstamp_fname, st1, "rawdata","dat"], '.')
        #                 fopath = InputDict["dumppath"]*"/"*fname_out
        #                 open(fopath, "w") do io
        #                     for i in 1:length(Schan.x)
        #                         write(io, @sprintf("%12.8f\n", Schan.x[i]))
        #                     end
        #                 end
        #
        #                 #dump kurtosis trace
        #                 mkpath(InputDict["dumppath"])
        #                 tstamp_fname = replace(tstamp, ":" => "-")
        #                 fname_out = join([tstamp_fname, st1, "kurtosis","dat"], '.')
        #                 fopath = InputDict["dumppath"]*"/"*fname_out
        #                 open(fopath, "w") do io
        #                     for i in 1:length(S1.misc["kurtosis"])
        #                         write(io, @sprintf("%12.8f\n", S1.misc["kurtosis"][i]))
        #                     end
        #                 end
        #
        #                 #dump removed trace
        #                 tstamp_fname = replace(tstamp, ":" => "-")
        #                 fname_out = join([tstamp_fname, st1, "remdata","dat"], '.')
        #                 fopath = InputDict["dumppath"]*"/"*fname_out
        #                 open(fopath, "w") do io
        #                     for i in 1:length(S1.x)
        #                         write(io, @sprintf("%12.8f\n", S1.x[i]))
        #                     end
        #                 end
        #             end
        #
        #
        #             #it's not allowed to save this into binary;
        #             delete!(S1.misc, "eqtimewindow")
        #             delete!(S1.misc, "kurtosis")
        #
        #         else
        #             #download error found: save as it is.
        #             S1 = Schan
        #         end
        #
        #         SremEQ += S1
        #
        #     end
        #
        #
        #     # if some of SeisChannels in Stemp have a data, save temp file
        #     y, jd = parse.(Int64, split(InputDict["DLtimestamplist"][dlid], ".")[1:2])
        #     tstamp_fname = replace(tstamp, ":" => "-")
        #     fname_out = join([tstamp_fname, st1, "FDSNWS","dat"], '.')
        #     # save as intermediate binary file
        #     t_write = @elapsed SeisIO.wseis(InputDict["tmppath"]*"/"*fname_out, SremEQ)
        #
        # end

        # __precompile__()
        # module Remove_eq
        #
        # export detect_eq_kurtosis, stalta, remove_eq, s_whiten, s_whiten!
        #
        # using SeisIO, SeisNoise, FFTW, DSP, StatsBase, Plots, JLD2, SeisIO, Printf
        #
        # """
        #     detect_eq_kurtosis(data::SeisChannel,tw::Float64=60.0, threshold::Float64=3.0, overlap::Float64=30)
        #
        # find earthquake by kurtosis threshold
        #
        # # Input:
        #     - `data::SeisChannel`    : SeisData from SeisIO
        #     - `tw::Float64`  : time window to evaluate if earthquake is contained.
        #     - `threshold::Float64` : kurtosis threshold: if kurtosis > threshold, the time window contains earthquake
        #     - `overlap::Float64`           : overlap of time window to control the margin of earthquake removal. (large overlap assigns large margin before and after earthquake.)
        #
        #     kurtosis evaluation following Baillard et al.(2013)
        # """
        # function detect_eq_kurtosis(data::SeisChannel; tw::Float64=60.0, kurtosis_threshold::Float64=3.0, overlap::Float64=30)
        #
        #     #convert window lengths from seconds to samples
        #     twsize = trunc(Int, tw * data.fs)
        #     overlapsize = trunc(Int, overlap * data.fs)
        #
        #     #calculate how much to move beginning of window each iteration
        #     slide = twsize-overlapsize
        #
        #     #kurtosis of timeseries
        #     ku1 = data.misc["kurtosis"][:]
        #     #reset long window counter and triggers for current channel
        #     i = 0
        #
        #     #loop through current channel by sliding
        #     while i < length(ku1) - twsize
        #
        #         #check if last window and change long window length if so
        #         if length(ku1) - i < twsize
        #             twsize = length(ku1)-i
        #         end
        #
        #         #define chunk of data based on long window length and calculate long-term average
        #         twtrace = @views ku1[i+1:i+twsize]
        #
        #         if !isnothing(findfirst(x -> x > kurtosis_threshold, twtrace))
        #             #this time window includes earthquake
        #             for tt= i+1:i+twsize
        #                 data.misc["eqtimewindow"][tt] = false
        #             end
        #         end
        #
        #         #advance long window
        #         i = i + slide
        #
        #     end
        #
        #     return data
        #
        # end
        #
        #
        # """
        #     detect_eq_stalta(data::SeisChannel,tw::Float64=60.0, threshold::Float64=3.0, overlap::Float64=30)
        #
        # find earthquake and tremors by STA/LTA
        #
        # # Input:
        #     - `data::SeisChannel`  : SeisChannel from SeisIO
        #     - `longwin::Float64`   : long time window
        #     - `shortwin::Float64`  : short time window
        #     - `threshold::Float64` : STA/LTA threshold: if STA/LTA > threshold, the time window contains earthquake
        #     - `overlap::Float64`   : overlap of time window to control the margin of earthquake removal. (large overlap assigns large margin before and after earthquake.)
        #     - `stalta_absoluteclip::Float64`   : clip if maximum absolute value exceeds this number (for the purpose of removing incoherent noise)
        #
        #     original code written by Seth Olinger. For our purpose, overlap is applied for short time window
        # """
        # function detect_eq_stalta(data::SeisChannel,longWinLength::Float64, shortWinLength::Float64, threshold::Float64, overlap::Float64, stalta_absoluteclip::Float64,
        #                             InputDict::Dict, tstamp::String, st1::String; datagap_eps::Float64=1e-8)
        #
        #
        #     #convert window lengths from seconds to samples
        #     longWin = trunc(Int,longWinLength * data.fs)
        #     shortWin = trunc(Int,shortWinLength * data.fs)
        #     overlapWin = trunc(Int,overlap * data.fs)
        #     trace = @view data.x[:]
        #
        #     #buffer for sta/lta value
        #     staltatrace = zeros(length(trace))
        #
        #     #calculate how much to move beginning of window each iteration
        #     slide = longWin
        #
        #     #save weight
        #     eqweight = deepcopy(data.misc["eqtimewindow"])
        #     #print("before manipulation"); println(count(eqweight))
        #
        #     # manipulate weight to avoid underestimation of lta due to data gap
        #     # threshold zero signal by 0.001% of lta within the timeseries, which encompasses small computational error such as 1e-21
        #     ltaall = StatsBase.mean(abs.(trace))
        #     for i = 1:length(trace)
        #         if isless(abs(trace[i]), ltaall*datagap_eps)
        #             # this is regarded as zero signal (data gap)
        #             eqweight[i] = false
        #         end
        #     end
        #
        #     #print("after manipulation"); println(count(eqweight))
        #
        #     #define long window counter and empty trigger vector
        #     i = 1
        #
        #     #loop through current channel by sliding
        #     while i < length(trace)
        #
        #         #check if last window and change long window length if so
        #         if length(trace) - i < longWin
        #             longWin = length(trace)-i
        #         end
        #
        #         #define chunk of data based on long window length and calculate long-term average
        #         longTrace = trace[i:i+longWin]
        #
        #         lta = StatsBase.mean(abs.(longTrace), weights(eqweight[i:i+longWin]))
        #
        #         #reset short window counter
        #         n = 0
        #         breakflag = false
        #         #loop through long window in short windows
        #         while n <= longWin - shortWin
        #
        #             #define chunk of data based on short window length and calculate short-term average
        #             shortTrace = @view trace[i+n:i+n+shortWin]
        #             shortTraceWeight = @view eqweight[i+n:i+n+shortWin]
        #             #sta = mean(abs.(shortTrace))
        #             sta = StatsBase.mean(abs.(shortTrace), weights(shortTraceWeight))
        #             #sta = StatsBase.mean(abs.(shortTrace))
        #             stamax = maximum(abs.(shortTrace))
        #             #calculate sta/lta ration
        #             staLta = sta/lta
        #
        #             #record detection time if sta/lta ratio exceeds threshold
        #             if staLta > threshold || stamax > stalta_absoluteclip
        #                 #triggers[i+n] = 1
        #                 #this time window includes earthquake
        #                 for tt= i+n:i+n+shortWin
        #                     data.misc["eqtimewindow"][tt] = false
        #                 end
        #             end
        #
        #             if breakflag
        #                 break;
        #             end
        #
        #             #advance short window
        #             n = n + shortWin - overlapWin
        #
        #             #adjust for the last edge of long window
        #             if n >= longWin - shortWin
        #                 n = longWin - shortWin
        #                 breakflag = true
        #             end
        #
        #             staltatrace[i+n+shortWin] = staLta
        #
        #         end
        #
        #         #advance long window
        #         i = i + slide
        #
        #     end
        #
        #     if InputDict["dumptraces"]
        #         #dump sta/lta trace
        #         mkpath(InputDict["dumppath"])
        #         tstamp_fname = replace(tstamp, ":" => "-")
        #         fname_out = join([tstamp_fname, st1, "stalta","dat"], '.')
        #         fopath = InputDict["dumppath"]*"/"*fname_out
        #         open(fopath, "w") do io
        #             for i in 1:length(staltatrace)
        #                 write(io, @sprintf("%12.8f\n", staltatrace[i]))
        #             end
        #         end
        #     end
        #     return data
        #
        # end
        #
        #
        # """
        #     remove_eq(data::SeisChannel)
        #
        # remove earthquake by kurtosis and STA/LTA threshold
        #
        # # Input:
        #     - `data::SeisData`    : SeisData from SeisIO
        #
        # """
        # function remove_eq(data::SeisChannel, data_origin::SeisChannel, plot_kurtosis_α::Float64, max_taper_dur::Float64,
        #     plot_boxheight::Float64, plot_span::Int64, plot_fmt::String, fodir::String, tstamp::String, tvec::Array{Float64,1}, IsSaveFig::Bool)
        #
        #     eqidlist = data.misc["eqtimewindow"][:]
        #     nx = length(data.x)
        #
        #     i = 1
        #
        #     t1 = []
        #     t2 = []
        #     y1 = []
        #     y2 = []
        #
        #     tt1 = 0
        #     tt2 = 0
        #     tt3 = 0
        #     tt4 = 0
        #     tt5 = 0
        #
        #     while i <= nx
        #         if !eqidlist[i]
        #             push!(t1, tvec[i])
        #
        #             t1id = i
        #
        #             #find next id
        #             tt1 += @elapsed nexttrueid = findfirst(x -> x == true, eqidlist[t1id:end])
        #
        #             if isnothing(nexttrueid)
        #                 # all data is removed
        #                 t2id = nx
        #                 iinc = nx
        #             else
        #                 t2id = (t1id - 1) + nexttrueid - 1 # first false id to end false id
        #                 iinc = t2id - i + 1 # first true as end of false + 1
        #             end
        #
        #             push!(t2, tvec[t2id])
        #
        #             max_wintaper_duration = Int(data.fs*max_taper_dur)
        #
        #             # compute α given maximum taper length
        #             eq_length = t2id-t1id+1
        #             taper_length = 2*max_wintaper_duration
        #             tukey_length = eq_length + taper_length
        #             invert_tukey_α = taper_length/tukey_length
        #
        #             # define inverted tukey window
        #             tt2 += @elapsed invtukeywin = -tukey(Int(tukey_length), invert_tukey_α) .+ 1
        #
        #             # slice tukey window if it exceeds array bounds
        #             tt3 += @elapsed if t1id < max_wintaper_duration + 1
        #                 left_overflow = (max_wintaper_duration-t1id)+1
        #                 invtukeywin = @views invtukeywin[left_overflow+1:end]
        #                 # leftmost t
        #                 left = 1
        #             else
        #                 # full taper length
        #                 left = t1id - max_wintaper_duration
        #             end
        #
        #             tt4 += @elapsed if (nx - t2id + 1) < max_wintaper_duration + 1
        #                 right_overflow = (max_wintaper_duration-(nx - t2id + 1))+1
        #                 invtukeywin = @views invtukeywin[1:end-right_overflow]
        #                 # rightmost t
        #                 right = nx
        #             else
        #                 # full taper length
        #                 right = t2id + max_wintaper_duration
        #             end
        #
        #             # apply tukey window
        #             tt5 += @elapsed data.x[left:right] .*= invtukeywin
        #
        #             #boxsize
        #             push!(y1, -plot_boxheight)
        #             push!(y2, plot_boxheight)
        #
        #         else
        #             iinc = 1
        #         end
        #
        #         i += iinc
        #
        #     end
        #
        #     #println([tt1, tt2, tt3, tt4, tt5])
        #
        #     if IsSaveFig
        #
        #         figdir = joinpath(fodir, "fig_removalEQ")
        #         if !ispath(figdir); mkpath(figdir); end
        #
        #         Plots.plot(bg=:white)
        #         normalized_amp = 0.5 * maximum(data_origin.x[1:plot_span:end])
        #
        #         Plots.plot!(tvec[1:plot_span:end], data_origin.x[1:plot_span:end]./ normalized_amp,
        #                     label="raw data", color="black")
        #
        #         Plots.plot!(tvec[1:plot_span:end], data.x[1:plot_span:end]./ normalized_amp,
        #                     label="after removal", color="blue")
        #
        #         #to plot kurtosis computed points
        #         kurtid = findall(x -> !iszero(data.misc["kurtosis"][x]), 1:length(data.misc["kurtosis"]))
        #         if isempty(filter(!isnan, data.misc["kurtosis"][kurtid]))
        #             kurtnormalize = 1.0
        #         else
        #             kurtnormalize = maximum(filter(!isnan, data.misc["kurtosis"][kurtid]))
        #         end
        #
        #         Plots.plot!(tvec[kurtid], data.misc["kurtosis"][kurtid] ./ kurtnormalize .* plot_kurtosis_α,
        #                     color="red",
        #                     label="kurtosis")
        #
        #         p = Plots.plot!(xlabel = "[hours]",
        #                     ylabel = "[m/s]",
        #                     title =  @sprintf("%s %s", data.id, tstamp),
        #                     size = (1200, 800))
        #
        #         figname = @sprintf("%s/%s_%s.%s", figdir, data.id, tstamp, plot_fmt)
        #         Plots.savefig(p, figname)
        #
        #
        #         # save jld2 fig file
        #         # trace1 = scatter(;x=tvec[1:plot_span:end], y=data_origin.x[1:plot_span:end]./ normalized_amp,
        #         #  mode="lines", name="raw data", line_color="black")
        #         #
        #         # trace2 = scatter(;x=tvec[1:plot_span:end], y=data.x[1:plot_span:end]./ normalized_amp,
        #         #  mode="lines", name="after remove", line_color="blue")
        #         #
        #         # #to plot kurtosis computed points
        #         # kurtid = findall(x -> !iszero(data.misc["kurtosis"][x]), 1:length(data.misc["kurtosis"]))
        #         #
        #         # if isempty(filter(!isnan, data.misc["kurtosis"][kurtid]))
        #         #     kurtnormalize = 1.0
        #         # else
        #         #     kurtnormalize = maximum(filter(!isnan, data.misc["kurtosis"][kurtid]))
        #         # end
        #         #
        #         # trace3 = PlotlyJS.scatter(;x=tvec[kurtid], y=data.misc["kurtosis"][kurtid] ./ kurtnormalize .* plot_kurtosis_α,
        #         #  line_color="red", mode="lines", name="kurtosis")
        #         #
        #         # if !isempty(t1)
        #         #  shapes = PlotlyJS.rect(t1, t2, y1, y2; fillcolor="#ff99cc", opacity=0.3, line_width=0)
        #         #  layout = PlotlyJS.Layout(shapes=shapes, width=1200, height=600,
        #         #      xaxis=attr(title="Time [hour]"),
        #         #      yaxis=attr(title="Normalized velocity"),
        #         #      font =attr(size=13),
        #         #      showlegend=true,
        #         #      title = @sprintf("%s %s", data.id, tstamp))
        #         #
        #         #      p = PlotlyJS.plot([trace1; trace2; trace3],layout)
        #         #
        #         # else
        #         #  layout = PlotlyJS.Layout(width=1200, height=600,
        #         #      xaxis=attr(title="Time [hour]"),
        #         #      yaxis=attr(title="Normalized velocity"),
        #         #      font =attr(size=12),
        #         #      showlegend=true,
        #         #      title = @sprintf("%s %s", data.id, tstamp))
        #         #
        #         #  p = PlotlyJS.plot([trace1; trace2; trace3],layout)
        #         #
        #         # end
        #         #
        #         # figdir = joinpath(fodir, "fig")
        #         # mkpath(figdir)
        #         # figname = @sprintf("%s/%s_%s.%s", figdir, data.id, tstamp, plot_fmt)
        #         # ORCA.savefig(p, figname)
        #         # #display(p)
        #         # #println("press return for next plot...")
        #         # #readline()
        #     end
        #
        #     return data
        #
        # end
        #
        # """
        #     s_whiten!(data::SeisChannel)
        #
        # Apply spectral whitening to seischannel trace
        # # Input:
        #     - `data::SeisChannel`    : SeisChannel from SeisIO
        #     - `freqmin::Real`: Pass band low corner frequency.
        #     - `freqmax::Real`: Pass band high corner frequency.
        #     - `pad::Int`: Number of tapering points outside whitening band.
        # """
        # function s_whiten!(data::SeisChannel, freqmin::Float64, freqmax::Float64;pad::Int=50)
        #
        #     # compute fft of time series
        #     FFT = rfft(data.x,1)
        #
        #     # to use SeisNoise.whiten!() prepare (N, 2) array
        #     FFTtemp = Array{Complex{Float32}}(undef, length(FFT), 2)
        #     FFTtemp[:, 1] = FFT
        #
        #     SeisNoise.whiten!(FFTtemp,freqmin,freqmax,data.fs, data.t[end, 1], pad=pad)
        #
        #     data.x = irfft(FFTtemp[:,1],data.t[end, 1],1)
        #     return nothing
        #
        # end
        # s_whiten(data::SeisChannel, freqmin::Float64, freqmax::Float64;pad::Int=50) = (U = deepcopy(data);
        #     s_whiten!(U,freqmin,freqmax,pad=pad);
        #     return U)
        #
        # end


        #
        # """
        # convert_tmpfile(InputDict::OrderedDict)
        #
        # convert temporal file in "./seisdownload_tmp" to prescribed format.
        # """
        # function convert_tmpfile_seisremoveeq(InputDict::OrderedDict)
        #
        # 	paths_all   = SeisIO.ls(InputDict["tmpdir_rem"])
        # 	fopath 		= joinpath(InputDict["fodir"], "EQRemovedData.jld2")
        # 	fo 			= jldopen(fopath, "w")
        #
        # 	# this is advanced mode to apply in order to isolate components at same stations
        # 	if InputDict["IsIsolateComponents"]
        # 		# isolate components based on priority dictionary
        # 		paths = isolate_components(paths_all, InputDict)
        # 	else
        # 		paths = paths_all
        # 	end
        #
        # 	varnamelist     = []
        #
        #     for path in paths
        #
        # 		S1 = jldopen(path, "r") do fi
        # 			fi["S"]
        # 		end
        #
        # 		varname = splitdir(path)[2][1:end-5]
        # 		groupname = join(split(varname, ".")[1:2], ".")
        #
        # 		# select output format
        #
        # 		if isempty(filter(x -> x==varname, varnamelist))
        # 			push!(varnamelist, varname)
        # 			fo[joinpath("Waveforms",groupname,varname)] = S1
        # 		end
        #
        # 		# remove tmpfile
        # 		rm(path)
        #     end
        #
        #     JLD2.close(fo)
        #
        # 	rm(InputDict["tmpdir_rem"], recursive=true, force=true)
        #
        #     return nothing
        # end
        #
        # """
        # isolate_components(paths_all::AbstractArray, InputDict::Dict)
        #
        # Isolate components based on InputDict["priority_channles"].
        #
        # I. What is the purpose?
        # 	Some stations have multiple stations with different channels
        # 	(e.g. BP.LCCB..BP1 and BP.LCCB.40.SP1 during transition.)
        # 	This increases number of xcorr processes, while the phisicall
        # 	meaning is auto correlation using different stations; thus we duplicate
        # 	auto-correlation in this case (e.g. BP1-BP1 and SP1-BP1 at same place).
        # 	To avoid that, we implemented this function to isolate components with
        # 	each stations based on priority, so that we can perform long-term dv/v
        # 	analysis without such duplications.
        #
        # II. Process flow
        #
        # 	1. Try to find channels which has same network.station, and same component
        # 	but different channel (let them as SP1 and BP1 in the network "BP")
        #
        # 	2. If found, search priority with `haskey(InputDict["priority_channles"], "BP")`.
        # 	If not found, arbitrually pick up the first one (SP1).
        #
        # 	3. To check if the priority is assigned to second one (BP1), search the order of it in Dictionary;
        # 	(e.g. findfirst("BP1", InputDict["priority_channles"]["BP"])). if not,
        # 	 pick up the first one (SP1).
        #
        # 	4. If the priority is assigned to second one (BP1), search the priority for first one (SP1);
        # 	(e.g. findfirst("SP1", InputDict["priority_channles"]["BP"])).
        # 	If not found,  pick up the second one (BP1) as it is listed in
        # 	priority dictionary.
        #
        # 	5. Compare the priority between first and second one, and pick up the
        # 	earlier one as representive channel at this station.
        #
        # III. Potential issue
        #
        # 	We assume that the physical measurement at same station does not change so much
        # 	that the cross-correlation is not influenced by the replacement of channels.
        # 	If it is not satisfied, it causes diference in the cross-correlatino result.
        # 	Please check the consistency between channels at same station if you find
        # 	some discontinuous result before and after switch of channel.
        #
        # """
        # function isolate_components(paths_all::AbstractArray, InputDict::Dict)
        #
        # 	iso_list = []
        #
        # 	for path in paths_all
        # 		tmp = split(path, "/")[end]
        # 		# read meta data from file name
        # 		ftmpname = split(tmp, ".")
        #
        # 		if occursin("-", ftmpname[3])
        # 			# format would be y, jd, T00-00-00, sta, loc, cha
        # 			y, d, tmpT, net, sta, loc, cha = split(ftmpname, ".")
        # 			#iso_stationinfo = (join([y, d, net, sta, loc], "-"), cha)
        # 			# (time, net, sta for find the station, channel name and component)
        # 			push!(iso_list, (path, join([y, d, tmpT, net, sta], "-"), net, sta, cha[1:2], cha[3]))
        #
        # 		else
        # 			@warn "Format of tmp file is not y, jd, time, sta, loc, cha. Turn off IsIsolateComponents."
        # 			return paths
        # 		end
        # 	end
        #
        # 	@show iso_list
        #
        # 	isocomp_idlist = []
        #
        # 	for (ista, current_sta) in enumerate(iso_list)
        # 		for jsta = ista:length(iso_list)
        # 			compared_sta = iso_list[j]
        # 			if current_sta[2] != compared_sta[2] || current_sta[6] != compared_sta[6]
        # 				# there is no conflict in channel, so add the current one to isocomp list
        # 				if ista ∉ isocomp_idlist
        # 					push!(isocomp_idlist, ista)
        # 				end
        # 			else
        # 				# here current_sta[2] == compared_sta[2] && current_sta[6] == compared_sta[6]
        # 				# perform process 2. we currently take into account priority with
        # 				# all stations in certain network
        #
        # 				iso_net = compared_sta[3]
        # 				current_cha = current_sta[5]
        # 				compared_cha = compared_sta[5]
        #
        # 				if !haskey(InputDict["priority_channles"],iso_net)
        # 					# second compared station does not have priority. take the current one
        # 					if ista ∉ isocomp_idlist
        # 						push!(isocomp_idlist, ista)
        # 					end
        #
        # 				else
        # 					# this has priority list; perform process 3. e.g. compared_sta[5] = "SP"
        # 					priority_compared = findfirst(compared_cha, InputDict["priority_channles"][iso_net])
        # 					if isempty(priority_compared)
        # 						# this channel has no priority. take the current one
        # 						if ista ∉ isocomp_idlist
        # 							push!(isocomp_idlist, ista)
        # 						end
        #
        # 					else
        # 						# this has priority so that search priority for current one.
        # 						priority_current = findfirst(current_cha, InputDict["priority_channles"][iso_net])
        # 						if isempty(priority_current)
        # 							# current channel has no priority. take the compared one
        # 							if jsta ∉ isocomp_idlist
        # 								push!(isocomp_idlist, jsta)
        # 							end
        #
        # 							# filter out ista if it's in iscomp_list
        # 							filter!(x -> x != ista, isocomp_idlist)
        #
        # 						else
        # 							# both current and compared one has priority. compare the order, and
        # 							# add or replace ista in isocomp_idlist
        # 							if priority_current > priority_compared
        # 								# take current one
        # 								if ista ∉ isocomp_idlist
        # 									push!(isocomp_idlist, ista)
        # 								end
        # 							else
        # 								# take compared one
        # 								if jsta ∉ isocomp_idlist
        # 									push!(isocomp_idlist, jsta)
        # 								end
        # 								# filter out ista if it's in iscomp_list
        # 								filter!(x -> x != ista, isocomp_idlist)
        # 							end
        # 						end
        # 					end
        # 				end
        # 			end
        # 		end
        # 	end
        #
        # 	#DEBUG:
        # 	for id in isocomp_idlist
        # 		temppath = paths_all[id]
        # 		tmp = split(temppath, "/")[end]
        # 		println(tmp)
        # 	end
        #
        # 	return paths_all[isocomp_idlist]
        #
        # end
