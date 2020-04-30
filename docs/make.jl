using Documenter, SeisMonitoring

makedocs(;
    modules=[SeisMonitoring],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/kura-okubo/SeisMonitoring.jl/blob/{commit}{path}#L{line}",
    sitename="SeisMonitoring.jl",
    authors="Kurama Okubo",
    assets=String[],
)

deploydocs(;
    repo="github.com/kura-okubo/SeisMonitoring.jl",
)
