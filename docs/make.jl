using RepertoireMetrics
using Documenter

DocMeta.setdocmeta!(RepertoireMetrics, :DocTestSetup, :(using RepertoireMetrics); recursive=true)

makedocs(;
    modules=[RepertoireMetrics],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com>",
    sitename="RepertoireMetrics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mashu.github.io/RepertoireMetrics.jl",
        edit_link="main",
        assets=String[],
        mathengine=Documenter.MathJax3()
    ),
    pages=[
        "Home" => "index.md",
        "Quick Start" => "quickstart.md",
        "Metrics" => "metrics.md",
        "Composable Metrics" => "composable.md",
        "API Reference" => "api.md",
    ],
    remotes=nothing,
    warnonly=[:missing_docs, :docs_block],
)

deploydocs(;
    repo="github.com/mashu/RepertoireMetrics.jl",
    devbranch="main",
)
