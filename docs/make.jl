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
        "Quantum Machine Learning" => "qml.md",
        "Quantum Error Correction" => [
            "Overview and encoded algorithms" => "qec.md",
            "Noisy fault-tolerant gadgets" => "fault_tolerance.md",
            "A fault-tolerant logical iSWAP" => "ft_iswap.md",
            "Codes, distances, and magic-state factories" => "magic_factory.md"
        ],
        "OpenQASM" => "qasm.md",
        "Functions" => "functions.md"
    ]
)

deploydocs(
    repo = "github.com/aydindeger/BlueTangle.jl.git",
    devbranch = "main",
    branch = "gh-pages",
    versions = nothing
)
