export call_param_general!
function call_param_general!(g::GtkGridLeaf, rowcount::Int)
        ParamList = [
                "DownloadType",
                "Start time",
                "End time",
                "Output Directory"
                ]

        for key in ParamList
                g[2,rowcount]= GtkLabel(key, selectable=true)    # Cartesian coordinates, g[x,y]
                g[3,rowcount]=try
                         GtkEntry(name=key, text=InputDict[key])   # Cartesian coordinates, g[x,y]
                catch
                @warn("Input file is missing the parameter $(key).")
                end

                signal_connect(g[3,rowcount], "changed") do widget
                        # println(get_gtk_property(widget, :name, String))
                        InputDict[get_gtk_property(widget, :name, String)] = get_gtk_property(widget, :text, String)
                end
                rowcount += 1
        end
        return rowcount
end
