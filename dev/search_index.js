var documenterSearchIndex = {"docs":
[{"location":"concept/#Development-concept-of-SeisMonitoring.jl","page":"Development concept","title":"Development concept of SeisMonitoring.jl","text":"","category":"section"},{"location":"concept/#Overview","page":"Development concept","title":"Overview","text":"","category":"section"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"We have developed SeisMonitoring.jl to perform the TB scale ambient seismic noise processings. The aims of this package are the following:","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"Mass downloading the continuous seismic data, removing transient signals, cross-correlations, stacking and measurement of dv/v.\nOptions of the processing steps such as the band-pass filtering using the wavelet transform, selective stacking and the dv/v measurement of stretching and MWCS.\nParalleling the processes using multi-node in cluster.\nOptimizing the use of memory and the frequency of I/O.","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"The kernels of SeisMonitoring.jl to handle the seismic waveforms are SeisIO.jl and SeisNoise.jl.","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"Name Functions Reference\nSeisIO.jl Handle the seismic waveforms Jones, J. P., Okubo, K., Clements, T., and Denolle, M. A. Seisio: a fast, efficient geophysical data architecture for the julia language. Seismol. Res. Lett., 91(4):2368–2377, 2020, doi:10.1785/0220190295.\nSeisNoise.jl Tools for ambient seismic noise processings Clements, T. and Denolle, M. A. Seisnoise.jl: ambient seismic noise cross correlation on the cpu and gpu in julia. Seismol. Res. Lett., 92(1):517–527, 2020, doi:10.1785/0220200192.\n  ","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"We handle the seismic waveforms using the data structures such as SeisIO.SeisData and SeisNoise.CorrData, which contains the meta data including the history of applied processes on the data, such as the filtering and tapering, and the waveforms.","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"The functions of SeisMonitoring.jl have been developed in the separated packages such as SeisDownload.jl for the ease of maintenance of the packages and dependencies, which is suitable for the Julia. However, we decided to merge them into a single package mainly to allow for using all the packages in one line, using SeisMonitoring.","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"The dependencies are managed by Julia system, and are listed in Project.toml. We minimized the dependencies, in particular not using the python-related modules and MPI.jl.","category":"page"},{"location":"concept/#Process-parallelization","page":"Development concept","title":"Process parallelization","text":"","category":"section"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"Most of the processes are parallelized using Julia-native Distributed.pmap function as the tasks associated with the ambient seismic noise processing can be process parallelized, i.e. asynchronized parallelization by station pairs, days, and frequency bands without the communication across the workers. The pmap allows for the multi-node parallelization. You can perform it by","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"using Distributed\nNP = 480 # for 56 cores * 10 nodes. Some margins of cores to increase the RAM per core.\naddprocs(NP) # Add the workers\nprint(\"After adding the procs: NP=$(nprocs())\");\n@everywhere using SeisMonitoring # You need to redefine the packages in all the processors.","category":"page"},{"location":"concept/#Parallelization-in-computing-cross-correlation","page":"Development concept","title":"Parallelization in computing cross-correlation","text":"","category":"section"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"The most computationally expensive process is to compute the cross-correlations. It is not only due to the number of tasks; the duplication of the file I/O needs to be optimized.","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"note: Example of cross-correlation in parallel\nLet us do the correlations across the stations STA1, STA2, STA3 after applying the FFT on the waveforms. Then, the tasks are STA1-STA2, STA1-STA3, STA2-STA3.If you read the seismic data or FFT files from the disk in each task, you need to read the STA1 twice from the disk. This is redundant, which could cause the damage in the disk.\nThe frequency of file I/O is very important to secure the scratch system and disks. See the documentation by TACC https://docs.tacc.utexas.edu/tutorials/managingio/ to optimize it.\nSo, you want to store all the FFTs in the memory and distribute them to the tasks. However, if you send the FFT from master core hosting the tasks to the workers in the different node, it takes time to transfer the data. You can run this metric, but it is inefficient comparing to conduct the tasks within a node.\nTherefore, we first parallelized the time window e.g. every 2 years, and submitted the jobs separately, and each job executes the cross-correlations of all possible station pairs parallelized with the number of cores in the node. This metric is optimized in the file I/O such that the file of seismic data is accessed only once, and is distributed via the memory within the node.\nIt is not practical to store the data in a single large file even in the HDF5-based JLD2 format. The I/O from a single large file is not recommended when you conduct the process parallelization. Therefore, we separated the seismic data into a daily length, and also output the correlation functions into the small pieces of '.jld2' files, which is gathered in the post-processing.","category":"page"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"We monitored the usage of CPUs and memory using TACC Remora . You can find the case study of scaling of parallelization in SeisMonitoring_Paper/Appx/Scaling_Frontera.","category":"page"},{"location":"concept/#Compatible-to-the-other-application","page":"Development concept","title":"Compatible to the other application","text":"","category":"section"},{"location":"concept/","page":"Development concept","title":"Development concept","text":"The default file format and processing parameters are tuned to process the data for the High Resolution Seismic Network (HRSN) (doi:10.7932/HRSN); however, the work flow can be generalized for any kind of the data set. You can also reuse the internal functions to filter the data. We left the comments and some deprecated functions as the reference for future development.  Please customize the functions under the MIT License.","category":"page"},{"location":"functions/#Key-functions","page":"Functions","title":"Key functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Here we list key functions used in the SeisMonitoring.jl. You can find the use of the functions in the tutorial of SeisMonitoring_Example.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"SeisMonitoring.init_project\nSeisMonitoring.set_parameter\nSeisMonitoring.get_parameter\nSeisMonitoring.run_job\nSeisMonitoring.assemble_seisdata\nSeisMonitoring.assemble_corrdata\nSeisMonitoring.cc_medianmute!\nSeisMonitoring.compute_frequency_decomposition\nSeisMonitoring.const_slice_codawindow","category":"page"},{"location":"functions/#SeisMonitoring.init_project","page":"Functions","title":"SeisMonitoring.init_project","text":"init_project((\n    ;\n    project_name::String = \"project\",\n    project_inputdir::String = \"./\",\n    project_outputdir::String = \"./\",\n    gui::String = true;\n    force::Bool=false\n)\n\nInitiate project directory where data is output. You can output data to local machine, external HDD, scratch, etc.\n\nArguments\n\nproject_name::String        : project name used as directory name [default: \"project\"]\nproject_inputdir::String    : absolute/relative path to make new input project directory [default: \".\"]\nproject_outputdir::String   : absolute/relative path to make new output project directory [default: \".\"]\ngui::Bool                   : true if you want to use gui (you can use this to initiate case studies)\nforce::Bool=false           : true if you want to remove existing file and init project.\n\nInput and output directories can be separated for the use on local/HDD/cloud/scratch file system.\n\nExamples\n\nrun project on your local machine\n\ninit_project(project_name=\"project_test\",\n             project_inputdir=\".\",\n             project_outputdir=\".\")\n\nrun project with external drive\n\ninit_project(project_name=\"project_external_drive\",\n             project_inputdir=\".\",\n             project_outputdir=\"/path-to-directory-in-your-external-drive\")\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.set_parameter","page":"Functions","title":"SeisMonitoring.set_parameter","text":"set_parameter(inputfilename::String, key::String, value::String)\n\nSet parameter of input dictionary in the inputfile.\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.get_parameter","page":"Functions","title":"SeisMonitoring.get_parameter","text":"get_parameter(inputfilename::String, parametername::String)\n\nget parameter from input file.\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.run_job","page":"Functions","title":"SeisMonitoring.run_job","text":"run_job(inputfile::String=\"\";\nseisdownload::Bool=true,\nseisremoveeq::Bool=true,\nseisxcorrelation::Bool=true,\nseisstack::Bool=true,\n)\n\nrunning job in the project folder.\n\nArguments\n\ninputfile::String       : absolute/relative path to input file (e.g. \"./project/inputfile/input.jl\")\n\nOptions\n\n'seisdownload::Bool'      : run seisdownload if true [default:true]\n'seisremoveeq::Bool'      : run seisremoveeq if true [default:true]\n'seisxcorrelation::Bool'  : run seisxcorrelation if true [default:true]\n'seisstack::Bool'         : run seisstack if true [default:true]\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.assemble_seisdata","page":"Functions","title":"SeisMonitoring.assemble_seisdata","text":"assemble_seisdata(C::SeisChannel)\n\nAssemble seisdata from starttime to endtime (target time window), importing data from jld2 with SeisMonitoring.jl format.\n\nSeisNoise.phase_shift!(SeisChannel) is NOT applied during this process.\n\nArguments\n\nnetstachan::String  : net.sta.loc.chan to be assembled\nfileio::JLD2.JLDFile: JLD2.JLDFile io of input jld2 file\nstarttime::DateTime : starttime to be assembled\nendtime::DateTime   : endtime to be assembled\ndata_contents_fraction::Float64=0.8 : if data exists more than data-contents-fraction within target window, return S.\n\nReturn\n\nS::SeisChannel: SeisChannel which contains data from starttime to endtime\n\nNote\n\nThis assemble function aims to realize robust input formats in different time scale and data situation.\n\nWarning: When marging the original data, tapering is applied to make the data continuous. This might cause an issue when The original data chunk is too small comparing with target window.\n\nBad example: save the data chunk every half an hour, and assemble the data into one day.\n\nTo avoid that, please set close enough between cc_time_unit (unit of target window length) and data chunk;\n\nGood example: save the data chunk every day, and assemble the data into one day or harf a day.\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.assemble_corrdata","page":"Functions","title":"SeisMonitoring.assemble_corrdata","text":"assemble_corrdata(C::SeisChannel)\n\nAssemble corrdata from starttime to endtime (target time window) within a frequency band, importing data from jld2 with SeisMonitoring.jl format.\n\nArguments\n\nfileio::JLD2.JLDFile: JLD2.JLDFile io of input jld2 file\n\n- stationpair::String : NOTE: This is deprecated because of the change of corrdata jld2 format. stationpair to be assembled\n\nstarttime::DateTime : starttime to be assembled\nendtime::DateTime   : endtime to be assembled\nfreqkey::String : Frequency band key used to decompose frequency contents of cc.\nmin_cc_datafraction::Float64 : minimum data fraction of cc within the request time.\nCorrData_Buffer::Dict : Dictionary of CorrData to optimize File IO.\nMAX_MEM_USE::AbstractFloat=4.0 : Maximum memory use; throw warning if the memory use exceeds this number.\nrename::Bool=true : true if rename C.name to nochan_stationpair.\n\nReturn\n\nC::CorrData: CorrData which contains data from starttime to endtime\n\nNote: 2021.06.07 adding \"rename::Bool=true\" to assemble corrdate with channel change. Rename to nochan_stationpair.\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.cc_medianmute!","page":"Functions","title":"SeisMonitoring.cc_medianmute!","text":"cc_medianmute(C::CorrData, cc_medianmute_max::Float64 = 10.0, cc_medianmute_min::Float64 = 0.0)\n\nMute cross-correlation function whose maximum amplitude is more than cc_medianmute_max*median(maximum(abs.(cross-correlation functions)))  and less than cc_medianmute_min*median((maximum(abs.(cross-correlation functions)))\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.compute_frequency_decomposition","page":"Functions","title":"SeisMonitoring.compute_frequency_decomposition","text":"compute_frequency_decomposition(C::CorrData, freqency_band::Array{Float64,1}; cc_bpfilt_method::String=\"Butterworth\",\n                             dj::Float64 = 1 / 12, α0::AbstractFloat = 0.0, αmax::AbstractFloat = 0.25)\n\ncompute frequency decomposition of cross-correlation function in corrdata\n\nArgument\n\nC::CoreData : CorrData contains broadband frequency contents.\ncc_bpfilt_method::String : Bandpass filtering method. \"Butterworth\" or \"Wavelet\"\ndj::AbstractFloat: Spacing between discrete scales. Default value is 1/12.\nα0::AbstractFloat=0.0: Lowest tapering fraction for frequency adaptive tapering.\nαmax::AbstractFloat=0.25: Highest tapering fraction for frequency adaptive tapering.\n\nReturn\n\nC_all::Array{CorrData, 1}: CorrData contains narrow frequency band\nfreqband::Array{Array{Float64,1},1} : frequency band used in the C_all.\n\n\n\n\n\n","category":"function"},{"location":"functions/#SeisMonitoring.const_slice_codawindow","page":"Functions","title":"SeisMonitoring.const_slice_codawindow","text":"const_slice_codawindow!(A, maxlag, fm, fs, dist, background_vel,\nmin_ballistic_twin, max_coda_length; zeropad=false)\n\nSlicing coda window based on prescribed velocity, maximum coda length and maximum coda length.\n\nAuthor: Kurama Okubo (https://github.com/kura-okubo) 2021.03.15\n\n\n\n\n\nconst_slice_codawindow!(A, maxlag, fm, fs, dist, background_vel,\nmin_ballistic_twin, max_coda_length; zeropad=false)\n\nSlicing coda window based on prescribed velocity, maximum coda length and maximum coda length.\n\nAuthor: Kurama Okubo (https://github.com/kura-okubo) 2021.03.15\n\n\n\n\n\n\n\n","category":"function"},{"location":"process_parameters/#Processing-parameters","page":"Processing parameters","title":"Processing parameters","text":"","category":"section"},{"location":"process_parameters/","page":"Processing parameters","title":"Processing parameters","text":"Here is the list of processing parameters defined in the mainparam.jl. You can configure the parameters, which is passed to SeisMonitoring.run_job ","category":"page"},{"location":"process_parameters/","page":"Processing parameters","title":"Processing parameters","text":"key default value type description\nproject_name project String project name.\nproject_inputdir ./project_INPUT String project input directory which you initiated with init_project().\nproject_outputdir ./project_OUTPUT String project output directory which you initiated with init_project().\nstarttime 2004-04-01T00:00:00 DateTime process start time\nendtime 2004-04-03T00:00:00 DateTime process end time\nsampling_frequency 20.0 Float64 [Hz] Processing sampling frequency. Downsampling is applied when downloaded original data is higher sampling frequency.\nfreqency_band [0.01, 0.1, 0.2, 0.5, 1.0, 2.0] Float64 Frequency bands to be analyzed.\nMAX_MEM_USE 3.0 Float64 [GB] Maximum memory use per core in the environment.\ndownload_time_unit 86400 Int64 [s] Unit time of data request. (e.g. request data for each 10minutes.\ndownload_margin 300 Int64 [s] Download margin to be clipped to avoid edge effect.\nrequeststation_file ./projectOUTPUT/defaultrequeststations.jld2 String Request station dataframe saved in JLD2. See default_requeststations.jld2\nIsResponseRemove true Bool True if remove instrumental response while downloading data.\nIsLocationBox false Bool True if using lat-lon box for request.\nreg [35.7, 36.1, -120.7, -120.2] Float64 minlat, maxlat, minlon, maxlon\nIstmpfilepreserved false Bool True if you want to preserve temporal files (same size as raw data.)\nIsXMLfilepreserved false Bool True if you want to preserve station xml files.\nnumstationperrequest 1 Int64 Advanced: number of station per one HTTP request.\noutputformat JLD2 String JLD2 or ASDF: use JLD2 if you perform the following processes with SeisMonitoring.jl\nRawData_path default String \"default\" or absolute/relative path to rawdata. \"default\" links to project output directory.\nIsKurtosisRemoval true Bool Apply Kurtosis removal.\nIsSTALTARemoval true Bool Apply STA/LTA removal.\nIsWhitening false Bool Apply Spectral whitening.\nfreqmin_whiten 0.1 Float64 Minimum cutoff frequency for spectral whitening\nfreqmax_whiten 1.0 Float64 Maximum cutoff frequency for spectral whitening\nAppend_alltraces false Bool Append kurtosis and stalta traces to SeisChannel (this increases data size)\nshorttime_window 180.0 Float64 Short-time window used to compute kurtosis and sta/lta\nlongtime_window 86400.0 Float64 Long-time window used to compute sta/lta\ntimewindow_overlap 60.0 Float64 Short-time window overlap to compute kurtosis and sta/lta\nkurtosis_threshold 2.0 Float64 Kurtosis removal threshold (The normal distribution of kurtosis is normalized to be zero.)\nstalta_threshold 1.2 Float64 STA/LTA removal threshold (For our purpose, this threshold is smaller than ordinal detection.)\nstalta_absoluteclip 0.1 Float64 [unit-of-data] clip the signal above this value (basically for instrumental error.)\nfixed_tukey_margin 30.0 Float64 [s] Fixed turkey margin; duration of decay outside of zero padding\nIsIsolateComponents false Bool Advanced: isolating comonents at same station\ncc_time_unit 86400 Int64 [s] Unit time of cross-correlation window. e.g. 60*60*24 = 86400 indicates daily-cross correlation.\ncc_len 3600 Int64 [s] short-time window cross-correlation length\ncc_step 1800 Int64 [s] cross-correlation window step\nmaxlag 100.0 Float64 [s] Maximum time lag of cross-correlation.\ncc_RawData_path default String \"default\" or absolute/relative path to rawdata. \"default\" links to project OUTPUT/EQRemovedData.jld2.\ncc_normalization deconvolution String none, coherence or deconvolution.\ncorr_type CC String Type of correlation: CC (standard cross-correlation) or PCC (phase cross-correlation). See also doc in SeisNoise.jl\npairs_option SubString{String}[\"11\", \"22\", \"33\"] Vector{String} \"all\" or list of component pairs. e.g. XX, YY, ZZ\nchanpair_type SubString{String}[\"all\"] Vector{String} \"all\" or list of channel pair type. e.g. auto-achan, cross-achan, cross-xchan\ndata_contents_fraction 0.8 Float64 Advanced: discard cross-correlation if data fraction within cc_time_unit is less that this value.\nIsOnebit false Bool Apply One-bit normalization.\nsmoothing_windowlength 7 Int64 Advanced: number of points for boxcar smoothing window on coherence and deconvolution.\nwater_level 0.0001 Float64 Advanced: waterlevel [0.0 if not applied] on spectrum normalization with coherence and deconvolution method.\ncc_bpfilt_method ButterWorth String Frequency decomposition method. \"Butterworth\" or \"Wavelet\".\ncc_taper_α0 0.1 Float64 Advanced: Lowest tapering fraction for frequency adaptive tapering.\ncc_taper_αmax 0.25 Float64 Advanced: Highest tapering fraction for frequency adaptive tapering.\ncc_medianmute_max 5.0 Float64 Advanced: Threshold factor of median mute within cc_time_unit. NCF is removed if maximum(abs.(corr[:,i])) > cc_medianmute_max * median(maximum(abs.(corr)), dims=1)\ncc_medianmute_min 0.1 Float64 Advanced: Threshold factor of median mute within cc_time_unit. NCF is removed if maximum(abs.(corr[:,i])) < cc_medianmute_min * median(maximum(abs.(corr)), dims=1)\nIsPreStack true Bool Advanced: Pre-stacking corrdata within each cc_time_unit when assembling the corrdata for the sake of saving memory use.\ntimechunk_increment 1 Int64 Advanced: Number of time chunk increment for parallelization: large number is more efficient, but increase memory use.\nstack_RawData_dir default String \"default\" or absolute/relative path to cc directory. \"default\" links to project OUTPUT/cc.\nuse_local_tmpdir true Bool True if using local /tmp diretory. Please set true when running in cluster to avoid massive file I/O.\nstack_method linear String stacking method: linear, selective, robust, pws, robustpws are available\ncollect_stationpairs true Bool true if correct station pairs. Stacking without this process does not work.\ncompute_reference true Bool true if compute reference stack for longterm stack.\ncompute_shorttimestack true Bool true if compute shorttime stack for continuous monitoring.\nstack_pairs_option SubString{String}[\"11\", \"22\", \"33\"] Vector{String} \"all\" or list of component pairs. e.g. XX, YY, ZZ\naveragestack_factor 1 Int64 Integer factor of cc_time_unit for stacking duration. e.g. cc_time_unit = 1day and averagestack_factor=30 provides 30days moving window average.\naveragestack_step 1 Int64 Step of averagestack window.\nmin_cc_datafraction 0.5 Float64 Advanced: discard cross-correlation if data fraction within stacking period is less that this value.\nreference_starttime 2004-04-01T00:00:00 DateTime reference start time\nreference_endtime 2004-04-02T00:00:00 DateTime reference end time\ndist_threshold 1.0 Float64 Threshold of distance used for selective stacking.\ndistance_type CorrDist String Advanced: Distance type used in selective stacking. See https://github.com/JuliaStats/Distances.jl for available types.\nIsZeropadBeforeStack false Bool Zero padding outside of coda window using tukey window before stacking.\nbackground_vel 2000.0 Float64 [m/s] Approximation of background wave velocity, just used for coda slicing.\nmin_ballistic_twin 1.0 Float64 [s] Explicit ballistic time window to remove coherence around zero timelag. This is aimed to remove it mainly for auto-correlation.\nmax_coda_length 60.0 Float64 [s] Maximum coda window length.\nmwcc_threshold 0.5 Float64 mwcc slice coda threshold.\nmwcc_len_α 3.0 Float64 moving window size factor (size = (mwcc_len_α/fm)*fs [point]).\nmin_codalength_α 1.0 Float64 Threshold of minimum codawindow length: min_codalength = min_codalength_α*mwcc window length.\ncoda_init_factor 2.0 Float64 [s] Coda window starts from coda_init_factor*dist/vel.\ncoda_minlen_factor 5.0 Float64 [s] Minimumlength is determined by this factor * (1/fm, period of cc) * fs points.\ncodaslice_debugplot false Bool If plot debug figures for coda slicing.\nnondim_max_coda_length 30.0 Float64 Deprecated: nondimensional maximum coda window length\nnondim_codamaxlag 60.0 Float64 Deprecated: coda max lag where kinetic energy is evaluated.\ncoda_energy_threshold -1.0 Float64 Deprecated: Advanced: Threshold for attenuation decay.\nIsAlternateRefChannel true Bool Advanced: Allow for using alternative station channel for reference. (e.g. BP.LCCB..BP1-BP.MMNB..BP1 is used as reference for BP.LCCB..SP1-BP.MMNB..SP1)\nkeep_corrtrace false Bool Advanced: Keep corr trace in CorrData if true. (require the strage to save corrs.)\nmeasurement_method mwcs String Stretching method for measuring dv/v and dQ^{-1}. \"stretching\",\"mwcs\",\"wcc\",\"dtw\",\"dualstretching\"\nmwcs_window_length 6.0 Float64 [s] The moving window length\nmwcs_window_step 3.0 Float64 [s] The step to jump for the moving window.\nmwcs_smoothing_half_win 5 Int64 [Points] MWCS smoothing half windown length.\nmwcs_max_dt 1.0 Float64 [s] MWCS threshold on dt.\nstretch_debugplot false Bool If plot debug figures for streching.\ndvv_stretching_range 0.02 Float64 Advanced: dvv stretching trial range for dvv (+- abs(dvv_stretching_range)).\ndvv_stretching_Ntrial 201 Int64 Advanced: dvv stretching trial number for dvv.\ngeometricalspreading_α 0.5 Float64 Advanced: geometrical spreading coefficient to compute Qcinv.\ncomputedqq_smoothing_windowlength 10.0 Float64 [s] smoothing windown length to compute envelope for compute_dvvdqq.","category":"page"},{"location":"#SeisMonitoring.jl","page":"Index","title":"SeisMonitoring.jl","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"This is the documentation of SeisMonitoring.jl","category":"page"},{"location":"#Installation","page":"Index","title":"Installation","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Type the commands below in the Julia REPL:","category":"page"},{"location":"","page":"Index","title":"Index","text":"using Pkg; Pkg.update();\nPkg.add(PackageSpec(name=\"SeisIO\", version=\"1.2.1\"));\nPkg.add(PackageSpec(name=\"SeisNoise\", version=\"0.5.3\"));\nPkg.add(url=\"https://github.com/kura-okubo/SeisDvv.jl\");\nPkg.add(url=\"https://github.com/kura-okubo/SeisMonitoring.jl\");","category":"page"},{"location":"","page":"Index","title":"Index","text":"note: Note\nThe metric above is a stable way to install the packages as the SeisDvv.jl and SeisMonitoring.jl are not registered to the General registory, and it may cause the issue when using ]add https://github.com/kura-okubo/SeisMonitoring.jl in the Julia REPL.","category":"page"},{"location":"","page":"Index","title":"Index","text":"tip: Uninstall the package\nType ] rm SeisMonitoring in the Julia REPL to remove the SeisMonitoring.","category":"page"},{"location":"#Tutorial","page":"Index","title":"Tutorial","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"We created the notebook of the tutorial in a different github repository, SeisMonitoring_Example. You can find how to download the data, remove the transient signals, compute cross-correlations, stack the correlation functions and measure the dv/v.","category":"page"},{"location":"","page":"Index","title":"Index","text":"You can access to the from the badge below:","category":"page"},{"location":"","page":"Index","title":"Index","text":"<a href=\"https://nbviewer.org/github/kura-okubo/SeisMonitoring_Example/blob/main/code/run_seismonitoring.ipynb\" target=\"_blank\">\n   <img align=\"left\"\n      src=\"https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.png\"\n      width=\"109\" height=\"20\">\n</a>\n<br><br>","category":"page"},{"location":"","page":"Index","title":"Index","text":"or read the QR code:","category":"page"},{"location":"","page":"Index","title":"Index","text":"<img src=\"./assets/QRcode_seismonitoring_example.png\" alt=\"QR\" width=\"150\"/>","category":"page"},{"location":"","page":"Index","title":"Index","text":"See also the repository of SeisMonitoring_Paper for the post-processing using the dv/v over 20 years at Parkfield.","category":"page"},{"location":"#Reference","page":"Index","title":"Reference","text":"","category":"section"}]
}
