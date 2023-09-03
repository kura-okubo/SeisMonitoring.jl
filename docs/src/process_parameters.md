# Processing parameters

Here is the list of processing parameters defined in the `mainparam.jl`. You can configure the parameters, which is passed to [`SeisMonitoring.run_job
`](@ref) 

| key | default value | type | description |
| :--- | :--- | :--- | :--- |
| project\_name | project | String | project name. |
| project\_inputdir | ./project_INPUT | String | project input directory which you initiated with init\_project(). |
| project\_outputdir | ./project_OUTPUT | String | project output directory which you initiated with init\_project(). |
| starttime | 2004-04-01T00:00:00 | DateTime | process start time |
| endtime | 2004-04-03T00:00:00 | DateTime | process end time |
| sampling\_frequency | 20.0 | Float64 | [Hz] Processing sampling frequency. Downsampling is applied when downloaded original data is higher sampling frequency. |
| freqency\_band | [0.01, 0.1, 0.2, 0.5, 1.0, 2.0] | Float64 | Frequency bands to be analyzed. |
| MAX\_MEM\_USE | 3.0 | Float64 | [GB] Maximum memory use per core in the environment. |
| download\_time\_unit | 86400 | Int64 | [s] Unit time of data request. (e.g. request data for each 10minutes. |
| download\_margin | 300 | Int64 | [s] Download margin to be clipped to avoid edge effect. |
| requeststation\_file | ./project_OUTPUT/default_requeststations.jld2 | String | Request station dataframe saved in JLD2. See default\_requeststations.jld2 |
| IsResponseRemove | true | Bool | True if remove instrumental response while downloading data. |
| IsLocationBox | false | Bool | True if using lat-lon box for request. |
| reg | [35.7, 36.1, -120.7, -120.2] | Float64 | minlat, maxlat, minlon, maxlon |
| Istmpfilepreserved | false | Bool | True if you want to preserve temporal files (same size as raw data.) |
| IsXMLfilepreserved | false | Bool | True if you want to preserve station xml files. |
| numstationperrequest | 1 | Int64 | Advanced: number of station per one HTTP request. |
| outputformat | JLD2 | String | JLD2 or ASDF: use JLD2 if you perform the following processes with SeisMonitoring.jl |
| RawData\_path | default | String | "default" or absolute/relative path to rawdata. "default" links to project output directory. |
| IsKurtosisRemoval | true | Bool | Apply Kurtosis removal. |
| IsSTALTARemoval | true | Bool | Apply STA/LTA removal. |
| IsWhitening | false | Bool | Apply Spectral whitening. |
| freqmin\_whiten | 0.1 | Float64 | Minimum cutoff frequency for spectral whitening |
| freqmax\_whiten | 1.0 | Float64 | Maximum cutoff frequency for spectral whitening |
| Append\_alltraces | false | Bool | Append kurtosis and stalta traces to SeisChannel (this increases data size) |
| shorttime\_window | 180.0 | Float64 | Short-time window used to compute kurtosis and sta/lta |
| longtime\_window | 86400.0 | Float64 | Long-time window used to compute sta/lta |
| timewindow\_overlap | 60.0 | Float64 | Short-time window overlap to compute kurtosis and sta/lta |
| kurtosis\_threshold | 2.0 | Float64 | Kurtosis removal threshold (The normal distribution of kurtosis is normalized to be zero.) |
| stalta\_threshold | 1.2 | Float64 | STA/LTA removal threshold (For our purpose, this threshold is smaller than ordinal detection.) |
| stalta\_absoluteclip | 0.1 | Float64 | [unit-of-data] clip the signal above this value (basically for instrumental error.) |
| fixed\_tukey\_margin | 30.0 | Float64 | [s] Fixed turkey margin; duration of decay outside of zero padding |
| IsIsolateComponents | false | Bool | Advanced: isolating comonents at same station |
| cc\_time\_unit | 86400 | Int64 | [s] Unit time of cross-correlation window. e.g. 60\*60\*24 = 86400 indicates daily-cross correlation. |
| cc\_len | 3600 | Int64 | [s] short-time window cross-correlation length |
| cc\_step | 1800 | Int64 | [s] cross-correlation window step |
| maxlag | 100.0 | Float64 | [s] Maximum time lag of cross-correlation. |
| cc\_RawData\_path | default | String | "default" or absolute/relative path to rawdata. "default" links to project OUTPUT/EQRemovedData.jld2. |
| cc\_normalization | deconvolution | String | none, coherence or deconvolution. |
| corr\_type | CC | String | Type of correlation: `CC` (standard cross-correlation) or `PCC` (phase cross-correlation). See also doc in SeisNoise.jl |
| pairs\_option | SubString{String}["11", "22", "33"] | Vector{String} | "all" or list of component pairs. e.g. XX, YY, ZZ |
| chanpair\_type | SubString{String}["all"] | Vector{String} | "all" or list of channel pair type. e.g. auto-achan, cross-achan, cross-xchan |
| data\_contents\_fraction | 0.8 | Float64 | Advanced: discard cross-correlation if data fraction within cc\_time\_unit is less that this value. |
| IsOnebit | false | Bool | Apply One-bit normalization. |
| smoothing\_windowlength | 7 | Int64 | Advanced: number of points for boxcar smoothing window on coherence and deconvolution. |
| water\_level | 0.0001 | Float64 | Advanced: waterlevel [0.0 if not applied] on spectrum normalization with coherence and deconvolution method. |
| cc\_bpfilt\_method | ButterWorth | String | Frequency decomposition method. "Butterworth" or "Wavelet". |
| cc\_taper\_α0 | 0.1 | Float64 | Advanced: Lowest tapering fraction for frequency adaptive tapering. |
| cc\_taper\_αmax | 0.25 | Float64 | Advanced: Highest tapering fraction for frequency adaptive tapering. |
| cc\_medianmute\_max | 5.0 | Float64 | Advanced: Threshold factor of median mute within cc\_time\_unit. NCF is removed if maximum(abs.(corr[:,i])) > cc\_medianmute\_max \* median(maximum(abs.(corr)), dims=1) |
| cc\_medianmute\_min | 0.1 | Float64 | Advanced: Threshold factor of median mute within cc\_time\_unit. NCF is removed if maximum(abs.(corr[:,i])) < cc\_medianmute\_min \* median(maximum(abs.(corr)), dims=1) |
| IsPreStack | true | Bool | Advanced: Pre-stacking corrdata within each cc\_time\_unit when assembling the corrdata for the sake of saving memory use. |
| timechunk\_increment | 1 | Int64 | Advanced: Number of time chunk increment for parallelization: large number is more efficient, but increase memory use. |
| stack\_RawData\_dir | default | String | "default" or absolute/relative path to cc directory. "default" links to project OUTPUT/cc. |
| use\_local\_tmpdir | true | Bool | True if using local /tmp diretory. Please set true when running in cluster to avoid massive file I/O. |
| stack\_method | linear | String | stacking method: linear, selective, robust, pws, robustpws are available |
| collect\_stationpairs | true | Bool | true if correct station pairs. Stacking without this process does not work. |
| compute\_reference | true | Bool | true if compute reference stack for longterm stack. |
| compute\_shorttimestack | true | Bool | true if compute shorttime stack for continuous monitoring. |
| stack\_pairs\_option | SubString{String}["11", "22", "33"] | Vector{String} | "all" or list of component pairs. e.g. XX, YY, ZZ |
| averagestack\_factor | 1 | Int64 | Integer factor of cc\_time\_unit for stacking duration. e.g. cc\_time\_unit = 1day and averagestack\_factor=30 provides 30days moving window average. |
| averagestack\_step | 1 | Int64 | Step of averagestack window. |
| min\_cc\_datafraction | 0.5 | Float64 | Advanced: discard cross-correlation if data fraction within stacking period is less that this value. |
| reference\_starttime | 2004-04-01T00:00:00 | DateTime | reference start time |
| reference\_endtime | 2004-04-02T00:00:00 | DateTime | reference end time |
| dist\_threshold | 1.0 | Float64 | Threshold of distance used for selective stacking. |
| distance\_type | CorrDist | String | Advanced: Distance type used in selective stacking. See https://github.com/JuliaStats/Distances.jl for available types. |
| IsZeropadBeforeStack | false | Bool | Zero padding outside of coda window using tukey window before stacking. |
| background\_vel | 2000.0 | Float64 | [m/s] Approximation of background wave velocity, just used for coda slicing. |
| min\_ballistic\_twin | 1.0 | Float64 | [s] Explicit ballistic time window to remove coherence around zero timelag. This is aimed to remove it mainly for auto-correlation. |
| max\_coda\_length | 60.0 | Float64 | [s] Maximum coda window length. |
| mwcc\_threshold | 0.5 | Float64 |  mwcc slice coda threshold. |
| mwcc\_len\_α | 3.0 | Float64 | moving window size factor (size = (mwcc\_len\_α/fm)\*fs [point]). |
| min\_codalength\_α | 1.0 | Float64 | Threshold of minimum codawindow length: min\_codalength = min\_codalength\_α\*mwcc window length. |
| coda\_init\_factor | 2.0 | Float64 | [s] Coda window starts from coda\_init\_factor\*dist/vel. |
| coda\_minlen\_factor | 5.0 | Float64 | [s] Minimumlength is determined by this factor \* (1/fm, period of cc) \* fs points. |
| codaslice\_debugplot | false | Bool | If plot debug figures for coda slicing. |
| nondim\_max\_coda\_length | 30.0 | Float64 | Deprecated: nondimensional maximum coda window length |
| nondim\_codamaxlag | 60.0 | Float64 | Deprecated: coda max lag where kinetic energy is evaluated. |
| coda\_energy\_threshold | -1.0 | Float64 | Deprecated: Advanced: Threshold for attenuation decay. |
| IsAlternateRefChannel | true | Bool | Advanced: Allow for using alternative station channel for reference. (e.g. BP.LCCB..BP1-BP.MMNB..BP1 is used as reference for BP.LCCB..SP1-BP.MMNB..SP1) |
| keep\_corrtrace | false | Bool | Advanced: Keep corr trace in CorrData if true. (require the strage to save corrs.) |
| measurement\_method | mwcs | String | Stretching method for measuring dv/v and dQ^{-1}. "stretching","mwcs","wcc","dtw","dualstretching"  |
| mwcs\_window\_length | 6.0 | Float64 | [s] The moving window length |
| mwcs\_window\_step | 3.0 | Float64 | [s] The step to jump for the moving window. |
| mwcs\_smoothing\_half\_win | 5 | Int64 | [Points] MWCS smoothing half windown length. |
| mwcs\_max\_dt | 1.0 | Float64 | [s] MWCS threshold on dt. |
| stretch\_debugplot | false | Bool | If plot debug figures for streching. |
| dvv\_stretching\_range | 0.02 | Float64 | Advanced: dvv stretching trial range for dvv (+- abs(dvv\_stretching\_range)). |
| dvv\_stretching\_Ntrial | 201 | Int64 | Advanced: dvv stretching trial number for dvv. |
| geometricalspreading\_α | 0.5 | Float64 | Advanced: geometrical spreading coefficient to compute Qcinv. |
| computedqq\_smoothing\_windowlength | 10.0 | Float64 | [s] smoothing windown length to compute envelope for compute\_dvvdqq. |
