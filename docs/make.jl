using Documenter
using HITRAN

makedocs(    
    modules = [HITRAN],
    authors="Oliver Kliebisch <oliver@kliebisch.net> and contributors",
    repo="https://github.com/tachawkes/hitran.jl/blob/{commit}{path}#L{line}",
    sitename="HITRAN.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Introduction & Motivation" => "index.md",
        "Quickstart/Example" => "quickstart.md",
        "Local HITRAN database" => "database.md",
        "Calculating spectra" => "spectra.md"        
    ],
)

deploydocs(
    repo = "github.com/TacHawkes/HITRAN.jl.git"
)
