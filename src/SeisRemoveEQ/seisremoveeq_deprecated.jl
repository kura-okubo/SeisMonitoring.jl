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
