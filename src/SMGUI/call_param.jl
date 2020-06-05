# This call_param_XXX!() is used to add parameters on GUI.
# Please add your new parameters into Defaultproject/set_default_inputdict.jl, then
# add the key name to ParamList.

function call_param_general!(g::GtkGridLeaf, rowcount::Int)
        ParamList = [
                "project_name",
                "project_inputdir",
                "project_outputdir",
                "starttime",
                "endtime",
                "sampling_frequency",
                "freqency_band",
                # "NP",
                "MAX_MEM_USE"
                ]

        for key in ParamList
                g[1,rowcount]= GtkLabel(key, selectable=true)    # Cartesian coordinates, g[x,y]
                g[2,rowcount]=try
                         GtkEntry(name=key, text=InputDict[key][1])   # Cartesian coordinates, g[x,y]
                catch
                        @warn("Your InputDict is missing the parameter $(key).")
                        rowcount += 1
                        continue
                end

                g[3,rowcount]= GtkLabel(string(InputDict[key][2]), selectable=true)    # Cartesian coordinates, g[x,y]
                g[4,rowcount]= GtkLabel(string(InputDict[key][3]), selectable=true)    # Cartesian coordinates, g[x,y]

                signal_connect(g[2,rowcount], "changed") do widget
                        # println(get_gtk_property(widget, :name, String))
                        replaced_key = get_gtk_property(widget, :name, String)
                        new_val = get_gtk_property(widget, :text, String)
                        set_values_inputdict!(InputDict, replaced_key, new_val)
                        # InputDict[get_gtk_property(widget, :name, String)] = get_gtk_property(widget, :text, String)
                end
                rowcount += 1
        end
        return rowcount
end

function call_param_seisdownload!(g::GtkGridLeaf, rowcount::Int)
        ParamList = [
                "download_time_unit",
                "download_margin",
                "requeststation_file",
                "IsResponseRemove",
                "IsLocationBox",
                "reg",
                "Istmpfilepreserved",
                "IsXMLfilepreserved",
                "outputformat"
                ]

        for key in ParamList
                g[1,rowcount]= GtkLabel(key, selectable=true)    # Cartesian coordinates, g[x,y]
                g[2,rowcount]=try
                         GtkEntry(name=key, text=InputDict[key][1])   # Cartesian coordinates, g[x,y]
                catch
                        @warn("Your InputDict is missing the parameter $(key).")
                        rowcount += 1
                        continue
                end

                g[3,rowcount]= GtkLabel(string(InputDict[key][2]), selectable=true)    # Cartesian coordinates, g[x,y]
                g[4,rowcount]= GtkLabel(string(InputDict[key][3]), selectable=true)    # Cartesian coordinates, g[x,y]

                signal_connect(g[2,rowcount], "changed") do widget
                        # println(get_gtk_property(widget, :name, String))
                        replaced_key = get_gtk_property(widget, :name, String)
                        new_val = get_gtk_property(widget, :text, String)
                        set_values_inputdict!(InputDict, replaced_key, new_val)
                        # InputDict[get_gtk_property(widget, :name, String)] = get_gtk_property(widget, :text, String)
                end
                rowcount += 1
        end
        return rowcount
end

function call_param_seisremoveeq!(g::GtkGridLeaf, rowcount::Int)
        ParamList = [
                "RawData_path",
                "IsKurtosisRemoval",
                "IsSTALTARemoval",
                "IsWhitening",
                "freqmin_whiten",
                "freqmax_whiten",
                "Append_alltraces",
                "shorttime_window",
                "longtime_window",
                "timewindow_overlap",
                "kurtosis_threshold",
                "stalta_threshold",
                "stalta_absoluteclip",
                "fixed_tukey_margin"
                ]

        for key in ParamList
                g[1,rowcount]= GtkLabel(key, selectable=true)    # Cartesian coordinates, g[x,y]
                g[2,rowcount]=try
                         GtkEntry(name=key, text=InputDict[key][1])   # Cartesian coordinates, g[x,y]
                catch
                        @warn("Your InputDict is missing the parameter $(key).")
                        rowcount += 1
                        continue
                end

                g[3,rowcount]= GtkLabel(string(InputDict[key][2]), selectable=true)    # Cartesian coordinates, g[x,y]
                g[4,rowcount]= GtkLabel(string(InputDict[key][3]), selectable=true)    # Cartesian coordinates, g[x,y]

                signal_connect(g[2,rowcount], "changed") do widget
                        # println(get_gtk_property(widget, :name, String))
                        replaced_key = get_gtk_property(widget, :name, String)
                        new_val = get_gtk_property(widget, :text, String)
                        set_values_inputdict!(InputDict, replaced_key, new_val)
                        # InputDict[get_gtk_property(widget, :name, String)] = get_gtk_property(widget, :text, String)
                end
                rowcount += 1
        end
        return rowcount
end

function call_param_seisxcorrelation!(g::GtkGridLeaf, rowcount::Int)
        ParamList = [
                "cc_time_unit",
                "cc_len",
                "cc_step",
                "maxlag",
                "cc_RawData_path",
                "cc_normalization",
                "corr_type",
                "pairs_option",
                "IsOnebit",
                "cc_bpfilt_method",
                "IsPreStack"
                ]

        for key in ParamList
                g[1,rowcount]= GtkLabel(key, selectable=true)    # Cartesian coordinates, g[x,y]
                g[2,rowcount]=try
                         GtkEntry(name=key, text=InputDict[key][1])   # Cartesian coordinates, g[x,y]
                catch
                        @warn("Your InputDict is missing the parameter $(key).")
                        rowcount += 1
                        continue
                end

                g[3,rowcount]= GtkLabel(string(InputDict[key][2]), selectable=true)    # Cartesian coordinates, g[x,y]
                g[4,rowcount]= GtkLabel(string(InputDict[key][3]), selectable=true)    # Cartesian coordinates, g[x,y]

                signal_connect(g[2,rowcount], "changed") do widget
                        # println(get_gtk_property(widget, :name, String))
                        replaced_key = get_gtk_property(widget, :name, String)
                        new_val = get_gtk_property(widget, :text, String)
                        set_values_inputdict!(InputDict, replaced_key, new_val)
                        # InputDict[get_gtk_property(widget, :name, String)] = get_gtk_property(widget, :text, String)
                end
                rowcount += 1
        end
        return rowcount
end

function call_param_seisstack!(g::GtkGridLeaf, rowcount::Int)
        ParamList = [
                "stack_RawData_dir",
                "stack_method",
                "compute_reference",
                "compute_shorttimestack",
                "stack_pairs_option",
                "averagestack_factor",
                "averagestack_step",
                "reference_starttime",
                "reference_endtime",
                "dist_threshold",
                "IsZeropadBeforeStack",
                "background_vel",
                "coda_Qinv",
                "min_ballistic_twin",
                "max_coda_length",
                # parameters for measurement
                "measurement_method",
                "stretch_distmethod"
                ]

        for key in ParamList
                g[1,rowcount]= GtkLabel(key, selectable=true)    # Cartesian coordinates, g[x,y]
                g[2,rowcount]=try
                         GtkEntry(name=key, text=InputDict[key][1])   # Cartesian coordinates, g[x,y]
                catch
                        @warn("Your InputDict is missing the parameter $(key).")
                        rowcount += 1
                        continue
                end

                g[3,rowcount]= GtkLabel(string(InputDict[key][2]), selectable=true)    # Cartesian coordinates, g[x,y]
                g[4,rowcount]= GtkLabel(string(InputDict[key][3]), selectable=true)    # Cartesian coordinates, g[x,y]

                signal_connect(g[2,rowcount], "changed") do widget
                        # println(get_gtk_property(widget, :name, String))
                        replaced_key = get_gtk_property(widget, :name, String)
                        new_val = get_gtk_property(widget, :text, String)
                        set_values_inputdict!(InputDict, replaced_key, new_val)
                        # InputDict[get_gtk_property(widget, :name, String)] = get_gtk_property(widget, :text, String)
                end
                rowcount += 1
        end
        return rowcount
end
