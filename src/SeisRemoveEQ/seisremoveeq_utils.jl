__precompile__()
module Utils

export convert_tmpfile, defaultinputdict!, printparams, initlogo

using SeisIO, Printf, Dates, JLD2, FileIO


"""
convert_tmpfile(InputDict::Dict)

convert temporal file in "./seisdownload_tmp" to prescribed format.
It has salvage mode, which allows to compile the temporal files in the case of failing during the download.
"""
function convert_tmpfile(InputDict::Dict; salvage::Bool=false)

    println("-------START CONVERTING-------")

	# save data to fopath file
	fopath = InputDict["fopath"]

	t = jldopen(InputDict["finame"])

	file = jldopen(fopath, "w")

	#!!!should be debuged because of Isocomponents!!!#
	file["info/stationlist"]     = t["info/stationlist"];

	if InputDict["IsStartendtime"]
		file["info/DLtimestamplist"] = InputDict["DLtimestamplist_selected"];
		file["info/starttime"]       = InputDict["starttime"];
		file["info/endtime"]         = InputDict["endtime"];
	else
		file["info/DLtimestamplist"] = t["info/DLtimestamplist"];
		file["info/starttime"]       = t["info/starttime"];
		file["info/endtime"]         = t["info/endtime"];
	end

	JLD2.close(t)

	# find all temporal files
    paths_all = ls(InputDict["tmppath"])
    fmt = InputDict["outputformat"]

    stationlist     = []
    DLtimestamplist = []
    varnamelist     = []

	if InputDict["IsIsolateComponents"]
		# isolate components based on priority dictionary
		paths = isolate_components(paths_all, InputDict)
	else
		paths = paths_all
	end

    for path in paths

        S = try
				SeisIO.rseis(path)[1]
			catch y
				println(y)
				@warn("cannot read tmpfile in seisremoveeq_tmp_sample. skipping")
				continue;
		end

        for ii = 1:S.n #loop at each seis channel

            # make station list
            staid = S[ii].id

            # save data (checking whether it's already in the jld2 because it causes an error)
            #parse info
            s_str = string(u2d(S[ii].t[1,2]*1e-6))

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

            else
                error("output format in $fmt is not implemented yet.")
            end
        end

		# remove tmpfile
		rm(path)

    end

    JLD2.close(file)

	rm(InputDict["tmppath"], recursive=true, force=true)

    return nothing
end

"""
isolate_components(paths_all::AbstractArray, InputDict::Dict)

Isolate components based on InputDict["priority_channles"].

I. What is the purpose?
	Some stations have multiple stations with different channels
	(e.g. BP.LCCB..BP1 and BP.LCCB.40.SP1 during transition.)
	This increases number of xcorr processes, while the phisicall
	meaning is auto correlation using different stations; thus we duplicate
	auto-correlation in this case (e.g. BP1-BP1 and SP1-BP1 at same place).
	To avoid that, we implemented this function to isolate components with
	each stations based on priority, so that we can perform long-term dv/v
	analysis without such duplications.

II. Process flow

	1. Try to find channels which has same network.station, and same component
	but different channel (let them as SP1 and BP1 in the network "BP")

	2. If found, search priority with `haskey(InputDict["priority_channles"], "BP")`.
	If not found, arbitrually pick up the first one (SP1).

	3. To check if the priority is assigned to second one (BP1), search the order of it in Dictionary;
	(e.g. findfirst("BP1", InputDict["priority_channles"]["BP"])). if not,
	 pick up the first one (SP1).

	4. If the priority is assigned to second one (BP1), search the priority for first one (SP1);
	(e.g. findfirst("SP1", InputDict["priority_channles"]["BP"])).
	If not found,  pick up the second one (BP1) as it is listed in
	priority dictionary.

	5. Compare the priority between first and second one, and pick up the
	earlier one as representive channel at this station.

III. Potential issue

	We assume that the physical measurement at same station does not change so much
	that the cross-correlation is not influenced by the replacement of channels.
	If it is not satisfied, it causes diference in the cross-correlatino result.
	Please check the consistency between channels at same station if you find
	some discontinuous result before and after switch of channel.

"""
function isolate_components(paths_all::AbstractArray, InputDict::Dict)

	iso_list = []

	for path in paths_all
		tmp = split(path, "/")[end]
		# read meta data from file name
		ftmpname = split(tmp, ".")

		if occursin("-", ftmpname[3])
			# format would be y, jd, T00-00-00, sta, loc, cha
			y, d, tmpT, net, sta, loc, cha = split(ftmpname, ".")
			#iso_stationinfo = (join([y, d, net, sta, loc], "-"), cha)
			# (time, net, sta for find the station, channel name and component)
			push!(iso_list, (path, join([y, d, tmpT, net, sta], "-"), net, sta, cha[1:2], cha[3]))

		else
			@warn "Format of tmp file is not y, jd, time, sta, loc, cha. Turn off IsIsolateComponents."
			return paths
		end
	end

	@show iso_list

	isocomp_idlist = []

	for (ista, current_sta) in enumerate(iso_list)
		for jsta = ista:length(iso_list)
			compared_sta = iso_list[j]
			if current_sta[2] != compared_sta[2] || current_sta[6] != compared_sta[6]
				# there is no conflict in channel, so add the current one to isocomp list
				if ista ∉ isocomp_idlist
					push!(isocomp_idlist, ista)
				end
			else
				# here current_sta[2] == compared_sta[2] && current_sta[6] == compared_sta[6]
				# perform process 2. we currently take into account priority with
				# all stations in certain network

				iso_net = compared_sta[3]
				current_cha = current_sta[5]
				compared_cha = compared_sta[5]

				if !haskey(InputDict["priority_channles"],iso_net)
					# second compared station does not have priority. take the current one
					if ista ∉ isocomp_idlist
						push!(isocomp_idlist, ista)
					end

				else
					# this has priority list; perform process 3. e.g. compared_sta[5] = "SP"
					priority_compared = findfirst(compared_cha, InputDict["priority_channles"][iso_net])
					if isempty(priority_compared)
						# this channel has no priority. take the current one
						if ista ∉ isocomp_idlist
							push!(isocomp_idlist, ista)
						end

					else
						# this has priority so that search priority for current one.
						priority_current = findfirst(current_cha, InputDict["priority_channles"][iso_net])
						if isempty(priority_current)
							# current channel has no priority. take the compared one
							if jsta ∉ isocomp_idlist
								push!(isocomp_idlist, jsta)
							end

							# filter out ista if it's in iscomp_list
							filter!(x -> x != ista, isocomp_idlist)

						else
							# both current and compared one has priority. compare the order, and
							# add or replace ista in isocomp_idlist
							if priority_current > priority_compared
								# take current one
								if ista ∉ isocomp_idlist
									push!(isocomp_idlist, ista)
								end
							else
								# take compared one
								if jsta ∉ isocomp_idlist
									push!(isocomp_idlist, jsta)
								end
								# filter out ista if it's in iscomp_list
								filter!(x -> x != ista, isocomp_idlist)
							end
						end
					end
				end
			end
		end
	end

	#DEBUG:
	for id in isocomp_idlist
		temppath = paths_all[id]
		tmp = split(temppath, "/")[end]
		println(tmp)
	end

	return paths_all[isocomp_idlist]

end
