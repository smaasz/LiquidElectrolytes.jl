## General remarks
All physical quantities are assumed to be consistently represented through their values expressed in basic SI units
(m, kg, s, A, K, mol, cd), supported by the [LessUnitful.jl](https://j-fu.github.io/LessUnitful.jl/) package
built on top of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).

## Electrolyte data

```@docs
ElectrolyteData
AbstractElectrolyteData
dlcap0(::ElectrolyteData)
debyelength(::ElectrolyteData)
chemical_potentials!
rrate
``` 

## Discretization system

```@docs
PNPSystem
pnpunknowns
electrolytedata
solventconcentration
```

## Standard calculations
```@docs
bulkbcondition
dlcapsweep
ivsweep
```

