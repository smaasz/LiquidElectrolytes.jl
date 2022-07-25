using Documenter, LiquidElectrolytes, LessUnitful


function mkdocs()
    DocMeta.setdocmeta!(LiquidElectrolytes, :DocTestSetup, :(using LiquidElectrolytes); recursive=true)
    makedocs(sitename="LiquidElectrolytes.jl",
             modules = [LiquidElectrolytes],
             clean = false, 
             doctest = true,
             authors = "J. Fuhrmann",
             repo="https://github.com/j-fu/LiquidElectrolytes.jl",
             pages=[
                 "Home"=>"index.md",
                 "API"=>"api.md",
                 "Internal API"=>"internal.md"
             ])
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/LiquidElectrolytes.jl.git", devbranch = "main")
    end

end

mkdocs()

