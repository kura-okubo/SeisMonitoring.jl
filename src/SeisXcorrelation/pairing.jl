export sort_pairs, get_corrtype
"""

    sort_pairs(pairs::AbstractArray)

Sort station pairs into their correlation category (auto-correlation, cross-correlation, cross-channel-correlation)

# Arguments
- `pairs::AbstractArray,`    : Array of station pairs, e.g. output of generate_pairs(file["info/stationlist"])

# Output
- `sorted_pairs::Dict{String, Array{String, N}}`    : dictionary containing station pairs grouped by correlation type

"""
function sort_pairs(pairs::AbstractArray)
    # dict containing each correlation category
    # preallocation of dict size assumes all stations have the same number of channels, which may not be true.
    sorted_pairs = Dict("auto-achan"     => [],
                                                "auto-xchan"     => [],
                                                "cross-achan"     => [],
                                                "cross-xchan"     => [])

    # fill dictionary based on detected correlation type
    for stnpair in pairs
        ct = get_corrtype(stnpair)
        push!(sorted_pairs[ct], stnpair)
    end
    #
    # # remove ["",""] used to initialize dict
    # sorted_pairs["acorr"]     = sorted_pairs["acorr"][:, 2:end]
    # sorted_pairs["xcorr"]     = sorted_pairs["xcorr"][:, 2:end]
    # sorted_pairs["xchancorr"] = sorted_pairs["xchancorr"][:, 2:end]

    return sorted_pairs
end

"""

    corrtype(stnpair::Array{String, 1})

Determine correlation into 4 types:
- 1. auto-achan: Same station, same channel e.g. X.A..HHZ-X.A..HHZ
- 2. auto-xchan: Same station, different channel e.g. X.A..HHZ-X.A..HHE
- 3. cross-achan: different station, same channel e.g. X.A..HHZ-X.B..HHZ
- 4. cross-xchan: different station, different channel e.g. X.A..HHZ-X.B..HHE

# Arguments
- `stnpair::Array{String, 1},`    : Station pair, e.g. ["BP.SMNB..BP1", "BP.SMNB..BP3"]

# Output
- `corrtype::String`    : correlation type, e.g. "xchancorr"

"""
function get_corrtype(stnpair::Array{String, 1})

    stn1 = join(split(stnpair[1], ".")[1:2], ".")
    cha1 = split(stnpair[1], ".")[4][end]
    stn2 = join(split(stnpair[2], ".")[1:2], ".")
    cha2 = split(stnpair[2], ".")[4][end]

    if (stn1 == stn2) && (cha1 == cha2)
        ct = "auto-achan"
    elseif (stn1 == stn2) && (cha1 != cha2)
        ct = "auto-xchan"
    elseif (stn1 != stn2) && (cha1 == cha2)
        ct = "cross-achan"
    elseif (stn1 != stn2) && (cha1 != cha2)
        ct = "cross-xchan"
    else
        warning("corrtype error with $(stn1).$(cha1)-$(stn2).$(cha2). at get_corrtype()")
    end

    return ct
#     if stnpair[1] == stnpair[2]
#         ct = "acorr"
#     # same station, different channel
# elseif (stnpair[1][end-3:end] != stnpair[2][end-3:end]) && (stnpair[1][1:end-3] == stnpair[2][1:end-3])
#         ct = "xchancorr"
#     # different station
#     else
#         ct = "xcorr"
#     end
    # return ct
end
