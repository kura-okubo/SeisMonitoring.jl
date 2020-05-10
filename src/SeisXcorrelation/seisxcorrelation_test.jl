using SeisIO, SeisNoise, JLD2, Dates, DataStructures
push!(LOAD_PATH,"/Users/kurama/.julia/dev/SeisMonitoring")
include("/Users/kurama/.julia/dev/SeisMonitoring/src/SeisXcorrelation/assemble_seisdata.jl")
include("/Users/kurama/.julia/dev/SeisMonitoring/src/SeisXcorrelation/seisxcorrelation_utils.jl")
fi = jldopen("/Volumes/Kurama_20190821/kurama/research/SeisMonitoring_dev/testproject_OUTPUT/seismicdata/EQRemovedData.jld2", "r")

starttime = DateTime("2004-04-01T00:00:00")
endtime = DateTime("2004-04-01T23:59:59.95")

sta1 = "BP.CCRB..BP1"
sta2 = "BP.EADB..BP1"
InputDict=Dict(
"data_contents_fraction" => 0.8,
"cc_len" => 3600,
"cc_step" => 1800,
"maxlag" => 120.0,
"freqency_band" => [0.1, 0.9],
"IsOnebit" => false,
"cc_method" => "cross-correlation",
"smoothing_half_win" => 5,
"waterlevel" => 1e-6,
"pairs_option" => lstrip.(split("all"))
)

StationDict = scan_stations("/Volumes/Kurama_20190821/kurama/research/SeisMonitoring_dev/testproject_OUTPUT/seismicdata/EQRemovedData.jld2")
# get station pairs
StationPairDict = get_stationpairs(StationDict, InputDict["cc_method"], InputDict["pairs_option"])


#1. assemble seisdata
t_assemble = @elapsed S1 , S2 = assemble_seisdata.([sta1, sta2], [fi], starttime, endtime, data_contents_fraction=InputDict["data_contents_fraction"])

#2. applying phase shift so that the data is aligned at the sampling frequency
phase_shift!.([S1, S2])

#3. convert to RawData
R1, R2 = RawData.([S1, S2], InputDict["cc_len"], InputDict["cc_step"])

#4. detrend, taper and band pass before computing fft and cross-correlation
clean_up!.([R1, R2], InputDict["freqency_band"][1], InputDict["freqency_band"][end])

#5. apply one-bit normalization if true
InputDict["IsOnebit"] && onebit!.([R1, R2])

#6. compute fft
t_fft = @elapsed FFT1, FFT2 = compute_fft.([R1, R2])

#7. spectral normalization
InputDict["cc_method"] == "coherence" && coherence!.([FFT1, FFT2], InputDict["smoothing_half_win"], InputDict["waterlevel"])

# first station is used as source; e.g. "BP.CCRB..BP1-BP.CCRB..BP1" then BP.CCRB..BP1 is used.
InputDict["cc_method"] == "deconvolution" && deconvolution!(FFT1, InputDict["smoothing_half_win"], InputDict["waterlevel"])

#8. Compute cross-correlation
t_xcorr = @elapsed xcorr = compute_cc(FFT1, FFT2, InputDict["maxlag"], corr_type=InputDict["cc_method"])
