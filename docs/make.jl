using Documenter, LiquidElectrolytes


#ENV["MPLBACKEND"]="agg"
#import PyPlot

function mkdocs()
    DocMeta.setdocmeta!(LiquidElectrolytes, :DocTestSetup, :(using LiquidElectrolytes); recursive=true)
    makedocs(sitename="LiquidElectrolytes.jl",
             modules = [LiquidElectrolytes],
             clean = false, 
             doctest = true,
             authors = "J. Fuhrmann",
             repo="https://github.com/j-fu/LiquidElectrolytes.jl",
             pages=[
                 "Home"=>"index.md"
             ])
end

mkdocs()

deploydocs(repo = "git@github.com:j-fu/LiquidElectrolytes.jl.git")

