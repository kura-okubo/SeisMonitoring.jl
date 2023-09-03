using Documenter, SeisMonitoring

makedocs(
	# format = Documenter.HTML(prettyurls = false),
    sitename = "SeisMonitoring.jl",
    authors = "kura-okubo",
    modules = [ SeisMonitoring ],
    pages = [
	"Index" => "index.md",
    "Development concept" => "concept.md",
	"Processing parameters" => "process_parameters.md",
	"Functions" => "functions.md"
    ])

deploydocs(
    repo = "github.com/kura-okubo/SeisMonitoring.jl.git",
)
