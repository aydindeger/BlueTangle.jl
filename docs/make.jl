push!(LOAD_PATH,"../src/")
import Pkg; Pkg.add("Documenter")
using Documenter
using BlueTangle

makedocs(
    sitename = "BlueTangle.jl",
    authors="Aydin Deger",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Advanced Examples" => "advanced.md",
        "Functions" => "functions.md"
    ]
)

deploydocs(
    repo = "github.com/aydindeger/BlueTangle.jl.git",
    devbranch = "main",
    branch = "gh-pages",
    versions = nothing
)