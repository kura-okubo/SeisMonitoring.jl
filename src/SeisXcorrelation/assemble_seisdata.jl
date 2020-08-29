const μs = 1.0e-6
const sμ = 1000000.0

"""
    assemble_seisdata(C::SeisChannel)

Assemble seisdata from starttime to endtime (target time window), importing data from
jld2 with SeisMonitoring.jl format.

SeisNoise.phase_shift!(SeisChannel) is **NOT** applied during this process.

# Arguments
- `netstachan::String`  : net.sta.loc.chan to be assembled
- `fileio::JLD2.JLDFile`: JLD2.JLDFile io of input jld2 file
- `starttime::DateTime` : starttime to be assembled
- `endtime::DateTime`   : endtime to be assembled
- `data_contents_fraction::Float64=0.8` : if data exists more than data-contents-fraction within target window, return S.

# Return
- `S::SeisChannel`: SeisChannel which contains data from starttime to endtime

# Note

This assemble function aims to realize robust input formats in different time scale and data situation.

**Warning**: When marging the original data, tapering is applied to make the data continuous.
This might cause an issue when **The original data chunk is too small comparing with target window.**

Bad example: save the data chunk every half an hour, and assemble the data into one day.

To avoid that, please set close enough between `cc_time_unit` (unit of target window length) and data chunk;

Good example: save the data chunk every day, and assemble the data into one day or harf a day.
"""
function assemble_seisdata(
    netstachan::String,
    fidir::String,
    starttime::DateTime,
    endtime::DateTime;
    data_contents_fraction::Float64 = 0.8,
)

    net, sta, loc, chan = split(netstachan, ".")
    #filtering file target
    files_all = readdir(fidir)
    files_station = filter(x -> occursin(netstachan, x), files_all)
    files_target = findall_target(files_station, starttime, endtime)

    if isempty(files_target)
        # the fileio does not have any file within starttime and endtime on the station
        return nothing
    end

    # read and merge seischannel
    S1 = SeisData(); removal_fraction_all = Float64[]
    bt_1 = @elapsed for file in files_target
        # Stemp = fileio[joinpath("Waveforms", join([net, sta], "."), file)]
        Stemp = rseis(joinpath(fidir, file))[1]
        haskey(Stemp.misc, "removal_fraction") && push!(removal_fraction_all, Stemp.misc["removal_fraction"])
        taper!(Stemp)
        S1 += Stemp
    end

    # the fileio does not have any file within starttime and endtime on the station
    isempty(S1) && return nothing

    # flatten data
    # S1 = merge(S1) # NOTE: to save memory use, replaced to merge!() to avoid deepcopy
    merge!(S1)
    # println(S1.t)
    # ungap at missing few sampling point due to rounding samplingcase
    # S1 = ungap(S1)
    ungap!(S1)
    # reconvert to seischannel

    #+++Deprecated the notes below; now using get_noisedatafraction to find noise data contents.+++#
    # NOTE: check the data amount per target window; if its fraction is less than data_contents_fraction,
    # discard it to avoid an issue with mean padding.
    # NOTE: If seisremoveeq is applied, the removal fraction is also taken into account.
    # datafrac_assemble = get_data_contents_fraction(S1[1],starttime,endtime)
    # datafrac_removeeq = (!isempty(removal_fraction_all) ? (1.0 - mean(removal_fraction_all)) : 1.0 )
    # datafraction_total = datafrac_assemble*datafrac_removeeq
    # println(datafraction_total)
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

    # syncrhonize seischannel so that it contains from starttime to endtime
    S1_sync = SeisIO.sync(S1, s=starttime, t=endtime, v=0)

    datafraction_total = get_noisedatafraction(S1_sync[1].x, zerosignal_minpts=100, eps_α=1e-6)

    if  datafraction_total < data_contents_fraction
        println("$(netstachan)-$(starttime)-$(endtime): data containts $(datafraction_total) is less than data_contents_fraction.")
        return nothing
    end
    # println("t_assemble: $(bt_1) at $(netstachan)-$(starttime)-$(endtime)")
    return S1_sync[1]
end

function findall_target(
    files_station::Array{String,1},
    starttime::DateTime,
    endtime::DateTime,
)

    files_target = []
    for file in files_station
        file_st, file_et = DateTime.(split(file, "__")[2:3])
        if !(starttime >= file_et || endtime <= file_st)
            # this file has overlap with target timewindow
            push!(files_target, file)
        end
    end
    return files_target
end

# function get_data_contents_fraction(
#     S1::SeisChannel,
#     starttime::DateTime,
#     endtime::DateTime,
# )
#     su, eu = SeisIO.t_win(S1.t, S1.fs) .* μs
#     target_window_length = (endtime - starttime).value / 1e3 # millisecond -> second
#     return (eu - su) / target_window_length
# end
