using Documenter
using TransitionPathTheory

makedocs(
    sitename = "TransitionPathTheory.jl",
    format = Documenter.HTML(),
    modules = [TransitionPathTheory],
    pages = [
        "Home" => "index.md",
        "Advanced Usage" => "advanced.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo = "github.com/70Gage70/TransitionPathTheory.jl.git",
    target = "build", 
    versions = nothing
)
