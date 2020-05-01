using Gtk
include("inputdict_io.jl")
include("default_param.jl")
include("param_general.jl")

"""
makeinput_gui()

Open Gtk gui to make input file for SeisMonitoring.jl
"""
function makeinput_gui(windowtitle::String="Default_param.jl")

  win = GtkWindow(windowtitle, 200, 40)

  #===Menu bar===#
  file = GtkMenuItem("_File")
  filemenu = GtkMenu(file)

  # open and load file
  open_ = GtkMenuItem("Open")
  push!(filemenu, open_)
  # save file
  save_ = GtkMenuItem("Save")
  push!(filemenu, save_)

  # quit gui
  push!(filemenu, GtkSeparatorMenuItem())
  quit_ = GtkMenuItem("Quit")
  push!(filemenu, quit_)

  mb = GtkMenuBar()
  push!(mb, file)  # notice this is the "File" item, not filemenu

  #===Parameter Entries===#
  add_row_overcol!(g::GtkGridLeaf, x::Int, e; colnum=1:4) = (g[colnum, x] = e; return x + 1)
  g = GtkGrid()
  rowcount = 1
  rowcount = add_row_overcol!(g, rowcount, mb)
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("General Setup"))
  rowcount = call_param_general!(g, rowcount)
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("SeisDownload"))
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("SeisRemoveEQ"))
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("SeisXcorrelation"))
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("SeisStack"))
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("SeisDvV"))
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("SeisDQ"))
  rowcount = add_row_overcol!(g, rowcount, GtkLabel("https://github.com/kura-okubo/SeisMonitoring.jl"))
  #

  set_gtk_property!(g, :column_homogeneous, true)
  set_gtk_property!(g, :row_spacing, 10)  # introduce a 15-pixel gap between columns
  set_gtk_property!(g, :column_spacing, 15)  # introduce a 15-pixel gap between columns
  push!(win, g)

  showall(win)

  # signal connection
  signal_connect(open_, "activate") do widget
    f = open_dialog("Open File", GtkNullContainer(), String[])
    # println(f)
    include(f)
    # println(InputDict)
    try
      makeinput_gui(f)
    catch
      @warn("Input file format error.")
    end
  end

  signal_connect(save_, "activate") do widget
    f = save_dialog("Save", GtkNullContainer(), String[])
    write_inputdict(f, InputDict)
    info_dialog("Input file is successfully saved.")
  end
  signal_connect(quit_, "activate") do widget
    ask_dialog("Do you want to quit window?", "No", "Yes") && Gtk.destroy(win)
  end

end
