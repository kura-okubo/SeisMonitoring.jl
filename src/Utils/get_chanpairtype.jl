"""

    get_chanpairtype(stnpair::Array{String, 1})

Determine correlation into 4 types:
- 1. auto-achan: Same station, same channel e.g. X.A..HHZ-X.A..HHZ
- 2. auto-xchan: Same station, different channel e.g. X.A..HHZ-X.A..HHE
- 3. cross-achan: different station, same channel e.g. X.A..HHZ-X.B..HHZ
- 4. cross-xchan: different station, different channel e.g. X.A..HHZ-X.B..HHE

# Arguments
- `stnpair::Array{String, 1},`    : Station pair, e.g. ["BP.SMNB..BP1", "BP.SMNB..BP3"]

# Output
- `corrtype::String`    : correlation type, e.g. "cross-xchan"

"""
function get_chanpairtype(stnpair::Array{String, 1})

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
        warning("corrtype error with $(stn1).$(cha1)-$(stn2).$(cha2). at get_chanpairtype()")
    end

    return ct

end
