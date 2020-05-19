# using Distributed
# addprocs(3)
#
using SeisIO, JLD2, DataStructures

#===The processes below are parallelized on each workers===#
#===(reading data and get group name)===#
function do_work(ch_paths, ch_seisdata) # define work function everywhere
	while true
	   path = take!(ch_paths)
	   # put!(ch_seisdata, tuple(splitdir(path)[2], myid()))
	   S = rseis(path)[1]
	   SC_all = [] # will contain all seischannels to be saved.

	   for ii = 1:S.n #loop at each seis channel
		   # make station list
		   staid = S[ii].id
		   # skip if S.t is empty
		   (isempty(S[ii]) || isempty(S[ii].t)) && continue;
		   s_str = string(u2d(S[ii].t[1,2]*1e-6))[1:19]
		   # compute end time:
		   et = S[ii].t[1,2]*1e-6 + (S[ii].t[end,1]-1)/S[ii].fs
		   e_str=string(u2d(et))[1:19]
		   # split network, station, location and channel
		   net, sta, loc, cha = split(S[ii].id, ".")
		   groupname = joinpath("Waveforms", join([net, sta], "."))
		   varname	  = join([net, sta, loc, cha], ".")*"__"*s_str*"__"*e_str*"__"*lowercase(cha)
		   push!(SC_all, (S[ii], groupname, varname))
	   end

	   put!(ch_seisdata, tuple(SC_all, myid()))
	end
end
#==========================================================================#

function convert_tmpfile(InputDict::OrderedDict, mode::String)

	paths_all = SeisIO.ls("./seisdownload_tmp")

	# this is advanced mode to apply in order to isolate components at same stations
	if InputDict["IsIsolateComponents"] && mode=="seisremoveeq"
		# isolate components based on priority dictionary
		paths = isolate_components(paths_all, InputDict)
	else
		paths = paths_all
	end


	ch_paths    = RemoteChannel(()->Channel{String}(length(paths)))
	ch_seisdata = RemoteChannel(()->Channel{Tuple}(Inf));
	n = length(paths);

	function make_jobs(n)
		for i in 1:n
			put!(ch_paths, paths[i])
		end
	end

	@async make_jobs(length(paths)); # feed the jobs channel with "n" jobs

	for p in workers() # start tasks on the workers to process requests in parallel
		remote_do(do_work, p, ch_paths, ch_seisdata)
	end

	# write SeisData into jld2 file at main process
	fodir   = InputDict["fodir"]
	fmt     = InputDict["outputformat"]
	if uppercase(fmt) == "JLD2"
		fopath = joinpath(fodir, "RawData.jld2")
		fo = jldopen(fopath, "w")
	elseif uppercase(fmt) == "ASDF"
		fopath = joinpath(fodir, "RawData.h5")
	else
		@error("outputformat: $(outputformat) is not available (JLD2 or ASDF).")
	end

	t_write = @elapsed while n > 0

		SC_all, where = take!(ch_seisdata) # take all seischannels from workers

		#===The data output process below should be done with single process.===#
		SCids = []
		for (SC, groupname, varname) in SC_all
			!haskey(fo, "Waveforms") && JLD2.Group(fo, "Waveforms")
			!haskey(fo, joinpath(groupname)) && JLD2.Group(fo, groupname)
			if !haskey(fo[groupname], varname)
				if lowercase(fmt) == "jld2"
					fo[joinpath(groupname,varname)] = SC
				elseif lowercase(fmt) == "asdf"
					write_hdf5(fopath, SC, add=true)
				end
			end
			push!(SCids, SC.id)
		end
		println("$(join(SCids, ":")) are processed on worker $where")
		#======================================================================#

		n = n - 1
	end

	uppercase(fmt) == "JLD2" && JLD2.close(fo)

	if !InputDict["Istmpfilepreserved"]
		rm("./seisdownload_tmp", recursive=true)
	end

	return t_write

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
