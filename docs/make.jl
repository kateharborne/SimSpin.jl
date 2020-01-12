push!(LOAD_PATH,"../src/")

using Documenter, SimSpin

makedocs(sitename="My Documentation")

deploydocs(
    repo = "github.com/kateharborne/SimSpin-Julia.git"
)
