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
                "NP"
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
                "savesamplefreq",
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
