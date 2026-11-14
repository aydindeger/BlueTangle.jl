using Documenter
using BlueTangle

makedocs(
    root = @__DIR__,
    sitename = "BlueTangle.jl",
    authors="Aydin Deger",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "QML" => "qml.md",
        "Quantum error correction" => "qec.md",
        "OpenQASM" => "qasm.md",
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
