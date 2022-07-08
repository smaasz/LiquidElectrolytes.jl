```@eval
using Markdown
Markdown.parse("""
$(read("../../README.md",String))
""")
```

## Data of Electrolyte

```@docs
ElectrolyteData
AbstractElectrolyteData
Cdl0(::ElectrolyteData)
ldebye(::ElectrolyteData)
chemical_potentials!
rrate
``` 

### Internal API
```@docs
LiquidElectrolytes.charge
LiquidElectrolytes.vrel
LiquidElectrolytes.c0_barc
LiquidElectrolytes.rlog
LiquidElectrolytes.wnorm
```

