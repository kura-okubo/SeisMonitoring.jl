module Utils
# with β version, please import SeisDownload.jl from the src directory as follows

include("downloadfunc.jl")
using .DownloadFunc

using SeisIO, Printf, Dates, JLD2, FileIO, Distributed

export get_starttimelist, get_timestamplist, get_stationlist, testdownload, convert_tmpfile, printparams, initlogo

"""

get_starttimelist(st::DateTime, et::DateTime, unittime::Float64)
calculate start time list for parallel downloading

    st: start time
    et: end time
    unittime: unit time in Second

    this function returns
    stlist: list of start time

    e.g.
    st = DateTime(2019,1,1,0,0,0)
    et = DateTime(2019,1,1,12,0,0)
    unittime = 3600

    stlist = get_starttimelist(st, et, unittime)

"""
function get_starttimelist(st::DateTime, et::DateTime, unittime::Real)

    reftime = st
    stlist = []

    while reftime < et
        push!(stlist, string(reftime))
        reftime += Dates.Second(float(unittime))
    end

    return stlist
end

"""
    get_timestamplist(stlist::Array{DateTime, 1})

    returns list of timestamp: Format = "Year_Julianday_Starttime".

"""
function get_timestamplist(stlist::Array{Any,1})

    timestamplist = []
    for stid = 1:length(stlist)
        yj = parse(Int64,stlist[stid][1:4])
        dj = md2j(yj, parse(Int64,stlist[stid][6:7]), parse(Int64,stlist[stid][9:10]))
        groupname    = string(yj)*"."*string(dj)*"."*stlist[stid][11:19] #Year_Julianday_Starttime
        push!(timestamplist, groupname)
    end

    return timestamplist

end

"""
    get_stationlist(network::Array{String, 1}, station::Array{String, 1}, location::Array{String, 1}, channel::Array{String, 1})

    returns list of request strings:

"""
function get_stationlist(network::Array{String, 1}, station::Array{String, 1}, location::Array{String, 1}, channel::Array{String, 1})

    stationlist = []
    for networkid = 1:length(network)
        for stationid = 1:length(station)
            for locationid = 1:length(location)
                for channelid = 1:length(channel)
                    requeststr = @sprintf("%s.%s.%s.%s", network[networkid], station[stationid], location[locationid], channel[channelid])
                    push!(stationlist, requeststr)
                end
            end
        end
    end

    return stationlist

end

"""
    testdownload(InputDict::Dict{String,Any} numofitr::Int64)

    print stats of download and return max_num_of_processes_per_parallelcycle

# Output
 -`max_num_of_processes_per_parallelcycle`: maximum number of processes for one request

"""
function testdownload(InputDict::Dict{String,Any}, numofitr::Int64)

    DownloadType    = InputDict["DownloadType"]

    trial_id          = 1
    InputDict_test = deepcopy(InputDict)
    test_suceededflag = false
    println("-------TEST DOWNLOAD START-----------")

    if DownloadType == "Noise" || DownloadType == "noise"

        while !test_suceededflag && trial_id < length(InputDict["starttimelist"])
            # select test request
            for j = 1:length(InputDict["stationinfo"]["stationlist"])
                InputDict_test["stationinfo"]["stationlist"] = [InputDict["stationinfo"]["stationlist"][j]]

                global t1 = @elapsed global dlerror = seisdownload_NOISE(trial_id, InputDict_test, testdownload=true) #[s]

                dl = [dlerror[i] for i in 1:length(dlerror)]
                if issubset(0, dl)
                    test_suceededflag = true
                    break;
                end
            end
            trial_id += 1
        end

    elseif  DownloadType == "Earthquake" || DownloadType == "earthquake"

        while !test_suceededflag && trial_id < length(InputDict["starttimelist"])
            # select test request
            for j = 1:length(InputDict["stationinfo"]["stationlist"])
                InputDict_test["stationinfo"]["stationlist"] = InputDict["stationinfo"]["stationlist"][j]

                global t1 = @elapsed global dlerror = seisdownload_EARTHQUAKE(trial_id, InputDict_test) #[s]

                dl = [dlerror[i] for i in 1:length(dlerror)]
                if issubset(0, dl)
                    test_suceededflag = true
                    break;
                end
            end
            trial_id += 1
        end
    end

    if !test_suceededflag
        error("All requests you submitted with input dictionary was failed. Please check the station availability in your request.")
    end

    estimated_downloadtime = now() + Second(round(3 * t1 * length(InputDict["DLtimestamplist"]) * numofitr / nprocs()))

    #println(mem_per_requestid)
    #println(max_num_of_processes_per_parallelcycle)
    println("-------DOWNLOAD STATS SUMMARY--------")

    println(@sprintf("Number of processes is %d.", nprocs()))

    # evaluate total download size by searching tmp directory
    tmppath = InputDict["tmppath"]
    s = read(`du -s -k $tmppath`, String)
    hdduse = parse(Int, split(s)[1])

    totaldownloadsize = hdduse * numofitr * length(InputDict["stationinfo"]["stationlist"])
    if totaldownloadsize < 1024 * 1024 # less than 1 GB
        totaldownloadsize = totaldownloadsize / 1024 #[MB]
        sizeunit = "MB"
    else
        totaldownloadsize = totaldownloadsize / 1024 / 1024
        sizeunit = "GB"
    end

    println(@sprintf("Total download size will be %4.2f [%s].", 0.8 * totaldownloadsize, sizeunit)) #0.8: considering compression efficiency
    println(@sprintf("Download will finish at %s.", round(estimated_downloadtime, Dates.Second(1))))
    println("*We have a time lag with downloading time above, like in 10 minutes or so.*")
    println("*This estimation also changes if some download requests fail and are skipped.*")
    println("-------START DOWNLOADING-------------")

    return nothing

end

"""
convert_tmpfile(InputDict::Dict)

convert temporal file in "./seisdownload_tmp" to prescribed format.
It has salvage mode, which allows to compile the temporal files in the case of failing during the download.
"""
function convert_tmpfile(InputDict::Dict; salvage::Bool=false)

    println("-------START CONVERTING-------------")
    paths = ls(InputDict["tmppath"])
    fopath = InputDict["fopath"]
    fmt = InputDict["outputformat"]


    #save HEADER information
    if fmt == "JLD2"
        file = jldopen(fopath, "w")
        file["info/DLtimestamplist"] = InputDict["DLtimestamplist"];
        file["info/starttime"]       = InputDict["starttime"]
        file["info/endtime"]         = InputDict["endtime"]
        file["info/DL_time_unit"]    = InputDict["DL_time_unit"]

    end

    stationlist     = []
    DLtimestamplist = []
    varnamelist     = []

    @simd for path in paths
        #println(path)
        try
            S = rseis(path)[1]
            #println(S)

            for ii = 1:S.n #loop at each seis channel

                # make station list
                staid = S[ii].id
                if isempty(filter(x -> x==staid, stationlist))
                    push!(stationlist, staid)
                end

                # save data (checking whether it's already in the jld2 because it causes an error)
                #parse info
                s_str = string(u2d(S[ii].t[1,2]*1e-6))

                # if time id is not matched with DLtimestamplist, this download is discarded for consistency even if it has some data.
                # here we also indirectly check if download_margin is correctly manipulated.
                if !isempty(filter(x -> x==s_str[1:19], InputDict["starttimelist"]))

                    # select output format
                    if fmt == "JLD2"
                        yj = parse(Int64, s_str[1:4])
                        mj = parse(Int64, s_str[6:7])
                        dj = parse(Int64, s_str[9:10])
                        tj = string(s_str)[11:19]

                        djm2j = md2j(yj, mj, dj)
                        groupname = string(yj)*"."*string(djm2j)*"."*tj #Year_Julianday_Starttime
                        varname = joinpath(groupname, staid)

                        if isempty(filter(x -> x==varname, varnamelist))
                            push!(varnamelist, varname)
                            file[varname] = S[ii]
                        end
                        # @info "save data $varname"

                    elseif fmt == "ASDF"
                        write_hdf5(fopath, S[ii])

                    else
                        error("output format in $fmt is not implemented yet.")
                    end

                end
            end

            if !InputDict["Istmpfilepreserved"]
    			rm(path)
    		end

        catch y
            #println(y)
        end
    end

    if fmt == "JLD2"
        # save final station list
        file["info/stationlist"] = stationlist
        JLD2.close(file)
    end


    return nothing
end


"""
printparams(param::Dict)

print parameters
"""
function printparams(param::Dict)
    printstyled("-----------Input Parameters-----------\n"; color=:cyan, bold=true)
    for key in keys(param)
        if length(string(param["$key"])) > 60
            param_str = string(param["$key"])[1:30]*"..."
        else
            param_str = string(param["$key"])
        end
        println(@sprintf("%-24s = %-10s", key, param_str))
    end
end


"""
initlogo()

print initial logo
"""
function initlogo()

    print("

      _____        _       _____                          _                    _
     / ____|      (_)     |  __ \\                        | |                  | |
    | (___    ___  _  ___ | |  | |  ___ __      __ _ __  | |  ___    __ _   __| |
     \\___ \\  / _ \\| |/ __|| |  | | / _ \\\\ \\ /\\ / /| '_ \\ | | / _ \\  / _` | / _` |
     ____) ||  __/| |\\__ \\| |__| || (_) |\\ V  V / | | | || || (_) || (_| || (_| |
    |_____/  \\___||_||___/|_____/  \\___/  \\_/\\_/  |_| |_||_| \\___/  \\__,_| \\__,_|
                      _         _  _
                     | |       | |(_)           |
    __      __       | | _   _ | | _   __ _     | v1.2 (Last update 08/15/2019)
    \\ \\ /\\ / /   _   | || | | || || | / _` |    | © Kurama Okubo
     \\ V  V /_  | |__| || |_| || || || (_| |    |
      \\_/\\_/(_)  \\____/  \\__,_||_||_| \\__,_|    |

")

    println("Job start running at "*string(now())*"\n")

end

end
