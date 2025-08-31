using Documenter
using MeshIntegrals

# Dynamically set all files in subdirectories of the source directory to include all files in these subdirectories.
# This way they don't need to be listed explicitly.
SPECIALIZATIONS_FILES = joinpath.(Ref("specializations"),
    readdir(joinpath(dirname(@__DIR__), "src", "specializations")))

makedocs(
    sitename = "MeshIntegrals.jl",
    pages = [
        "Home" => [
            "About" => "index.md",
            "Tutorial" => "tutorial.md",
            "Integration Rules" => "integration_rules.md",
            "Support Status" => "support.md",
            "Tips" => "tips.md"
        ],
        "Examples" => [
            "Changelog" => "examples/darts.md"
        ],
        "Developer Notes" => [
            "Changelog" => "developer/CHANGELOG.md",
            "How it Works" => "developer/how_it_works.md",
            "Specializations" => "developer/specializations.md"
        ],
        "Public API" => "api.md"
    ]
)

deploydocs(repo = "github.com/JuliaGeometry/MeshIntegrals.jl.git",
    devbranch = "main",
    push_preview = true)
