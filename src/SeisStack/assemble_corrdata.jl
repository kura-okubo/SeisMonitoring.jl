# include("get_cc_contents_fraction.jl")
using SeisIO, SeisNoise, Dates, JLD2, Statistics
"""
    assemble_corrdata(C::SeisChannel)

Assemble corrdata from starttime to endtime (target time window) within a frequency band, importing data from
jld2 with SeisMonitoring.jl format.

# Arguments
- `fileio::JLD2.JLDFile`: JLD2.JLDFile io of input jld2 file
# - `stationpair::String` : NOTE: This is deprecated because of the change of corrdata jld2 format. stationpair to be assembled
- `starttime::DateTime` : starttime to be assembled
- `endtime::DateTime`   : endtime to be assembled
- `freqkey::String` : Frequency band key used to decompose frequency contents of cc.
- `min_cc_datafraction::Float64` : minimum data fraction of cc within the request time.
- `CorrData_Buffer::Dict` : Dictionary of CorrData to optimize File IO.
- `MAX_MEM_USE::AbstractFloat=4.0` : Maximum memory use; throw warning if the memory use exceeds this number.
- `rename::Bool=true` : true if rename C.name to nochan_stationpair.

# Return
- `C::CorrData`: CorrData which contains data from starttime to endtime

Note: 2021.06.07 adding "rename::Bool=true" to assemble corrdate with channel change. Rename to nochan_stationpair.
"""
function assemble_corrdata(
    fileio,
    # stachanpair::String,    #BP.CCRB..BP1-BP.EADB..BP1
    starttime::DateTime,
    endtime::DateTime,
    freqkey::String;
    min_cc_datafraction::Float64 = 0.5,
    CorrData_Buffer::Dict=Dict(),
    MAX_MEM_USE::AbstractFloat=4.0, #[GB]
    rename::Bool=true # true if rename C.name to nochan_stationpair

)

    # corrdata_cputime debug
    dubug_t1 = 0; dubug_t2 = 0; dubug_t3 = 0; dubug_t4 = 0; dubug_t5 = 0;
    dubug_t6 = 0; dubug_t7 = 0;

    # Nfreqband = length(frequency_band) - 1
    # freqband = map(i -> [frequency_band[i], frequency_band[i+1]], 1:Nfreqband)
    C = CorrData()

    # 1. find all target timewindow
    cc_unit_time_all = keys(fileio)

    # find all path within target time window
    dubug_t1 = @elapsed files_target = findall_target_cc(cc_unit_time_all, starttime, endtime)
    isempty(files_target) && return C, CorrData_Buffer  # the fileio does not have any ccdata within starttime and endtime on the station. return empty C_all.

    # C_all = Dict{String, CorrData}()
    current_abskey_list = String[] # to be used to update CorrData_Buffer
    # 2. read and merge CorrData for all frequency band
    ccfracs = []
    for file in files_target # e.g. 2004-04-01T00:00:00--2004-04-02T00:00:00
        # abskey = joinpath(stachanpair, file, freqkey)
        abskey = joinpath(file, freqkey)

        if haskey(CorrData_Buffer, abskey)
            # read CorrData Buffer; which is prestacked when IsPreStack==true
            # println("debug: use buffer")
            dubug_t2 += @elapsed Ctemp = CorrData_Buffer[abskey]
            ccfrac = Ctemp.misc["ccfrac_within_cc_time_unit"]

        else
            # read data from file; performing Prestack if true and add to CorrData_Buffer
            # println("debug: read data")
            dubug_t3 += @elapsed Ctemp = try
                fileio[joinpath(abskey)]
            catch
                continue # skipping this time window because of no data
            end

            # (isnothing(Ctemp) || isempty(Ctemp)) && continue # NOTE: causing segmentation error at isempty(Ctemp)
            isempty(Ctemp.t) && continue

            if rename && (length(split(Ctemp.name, "-")) != 3)
                # rename to nochan_stationpair to avoid error in appending at C += Ctemp
                stachantemp = split(Ctemp.name, ".")
                # println(stachantemp)
                net1, sta1, loc1, cha1 = stachantemp[1:4]
                net2, sta2, loc2, cha2 = stachantemp[5:8]
                Ctemp.name = joinpath("$(net1).$(sta1)-$(net2).$(sta2)-$(cha1[end])$(cha2[end])")
            end

            # evaluate cc contents fraction
            # st, et = DateTime.(split(file, "--"))

            # NOTE: Ctemp.misc["ccfrac_within_cc_time_unit"] is computed at SeisXcorrelation stage.
            # dubug_t4 += @elapsed ccfrac = get_cc_contents_fraction(Ctemp,st,et)
            # Ctemp.misc["tccfrac_within_cc_time_unit"] = ccfrac
            ccfrac = Ctemp.misc["ccfrac_within_cc_time_unit"]

            # add CorrData to CorrData_Buffer after stacking.
            CorrData_Buffer[abskey] = Ctemp
            push!(current_abskey_list, abskey)
        end

        push!(ccfracs, ccfrac)

        C += Ctemp

        # check memory use
        memuse = round(sizeof(C.corr) * 1e-9, digits=6) #[GB]


        if memuse > MAX_MEM_USE
            @warn("Momory use will be $(memuse*Nfreqband) GB, which is more than MAX_MEM_USE=$(MAX_MEM_USE) GB.
                   This may cause memory overflow in your environment. Please increase MAX_MEM_USE.")
        end
    end
    # if C is nothing or empty, return empty CorrData
    # (isnothing(C) || isempty(C)) && (C = CorrData())
    isempty(C.t) && (C = CorrData())

    if  mean(ccfracs) < min_cc_datafraction
        # println("debug: data containts $(mean(ccfracs)) is less than cc_contents_fraction.")
        C = CorrData()
    end

    C.misc["cc_contents_frac"] = mean(ccfracs)

    # C_all[freqkey] = C1
    # end

    # update CorrData_Buffer
    dubug_t7 += @elapsed update_CorrData_Buffer!(CorrData_Buffer, current_abskey_list)

    return C, CorrData_Buffer
end

function findall_target_cc(
    cc_unit_time_all::Array{String,1},
    starttime::DateTime,
    endtime::DateTime,
)

    files_target = []
    for file in cc_unit_time_all
        file_st, file_et = DateTime.(split(file, "--"))
        if (starttime <= file_st) && (file_et <= endtime)
            # this file is within the request target window
            push!(files_target, file)
        end
    end
    return files_target
end

"""
    update_CorrData_Buffer!(CorrData_Buffer::Dict, current_abskey_list::Array{String, 1})

update CorrData_Buffer.

# Process flow

1. evaluate the symmetric difference of elements between abskey and current_abskey_list
2. if the symdiff_abskey is not in the current_abskey_list, delete that corrdata as it's no longer used in the process.
"""
function update_CorrData_Buffer!(CorrData_Buffer::Dict, current_abskey_list::Array{String, 1})
    symdiff_abskeys = symdiff(collect(keys(CorrData_Buffer)), current_abskey_list)

    for symdiff_abskey in symdiff_abskeys
        if symdiff_abskey ∉ current_abskey_list
            delete!(CorrData_Buffer, symdiff_abskey)
            # println("debug: $(symdiff_abskey) is removed from buffer.")
        end
    end
    return nothing
end

# Test
# using SeisIO, SeisNoise, JLD2, Dates
#
# fi = jldopen("/Volumes/Kurama_20190821/kurama/research/SeisMonitoring_dev/testproject_OUTPUT/cc/BP.CCRB-BP.EADB.jld2", "r")
# stachanpair = "BP.CCRB..BP1-BP.EADB..BP1"
# starttime = DateTime("2004-04-01T00:00:00")
# endtime = DateTime("2004-04-02T00:00:00")
# frequency_band = [0.01, 0.1, 0.2, 0.5, 1.0, 2.0]
# C, CorrData_Buffer = assemble_corrdata(fi,stachanpair,starttime,endtime,frequency_band, MAX_MEM_USE=2.0)
#
# println(CorrData_Buffer)
# # use corrdata buffer
# C, CorrData_Buffer = assemble_corrdata(fi,stachanpair,starttime,endtime,frequency_band,
#                         CorrData_Buffer=CorrData_Buffer, MAX_MEM_USE=0.005)
