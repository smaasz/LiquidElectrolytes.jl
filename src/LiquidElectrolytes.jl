module LiquidElectrolytes
using Base: @kwdef
using VoronoiFVM,ExtendableGrids
using DocStringExtensions
using CompositeStructs
using Parameters

# until the PR over there is accepted
CompositeStructs.to_expr(t::Number) = t


include("units.jl")

@phconstants N_A e R ε_0
const F=N_A*e
@siunits K mol dm m s g

export @siunits,@phconstants



include("electrolyte.jl")
export ElectrolyteData,Cdl0,c0_barc,chemical_potentials!, rrate

include("nppsystem.jl")
export NPPSystem
export nppunknowns,electrolytedata


include("cells.jl")
export voltagesweep,doublelayercap,bulkbc, TwoElectrodeCell



end # module
