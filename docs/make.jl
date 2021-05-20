push!(LOAD_PATH,"../src/")

using Documenter, SimSpin

makedocs(sitename="SimSpin.jl")

deploydocs(
    repo = "github.com/kateharborne/SimSpin.jl.git"
)
