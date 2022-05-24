````@eval
using Markdown
Markdown.parse("""
$(read("../../README.md",String))
""")
````

## Data of Electrolyte


```@docs
AbstractElectrolyteData
ElectrolyteData
``` 

### Internal API
```@docs
LiquidElectrolytes.charge
LiquidElectrolytes.vrel
```

## Physical unit handling

Physical unit handling is based on `Unitful.jl` and `PhysicalContstants.jl`,
prospectively this will be factorized out into another package. The rationale of these
macros is the fact that `Unitful.jl` facilitates working with strongly typed values.
Slightly counterintuitively to me, units with different prefixes are different types.
Moreover parts of the Julia infrastructure is not (yet) written in Julia, prominently
so sparse matrix solvers, and these are not able to handle units.

So we (mis)use the unit packages as databases and allow to retrive numerical values
in SI base units

```@docs
@si_str
@siunits
@phconstants
```
