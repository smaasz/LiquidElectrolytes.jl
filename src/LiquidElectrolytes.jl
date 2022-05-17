module LiquidElectrolytes
using VoronoiFVM,ExtendableGrids
using DocStringExtensions
using CompositeStructs
using Parameters
using ProgressLogging
using StaticArrays
using LinearAlgebra

include("units.jl")

@phconstants N_A e R ε_0
const F=N_A*e
@siunits K mol dm m s g

export @siunits,@phconstants



include("electrolyte.jl")
export ElectrolyteData,Cdl0,c0_barc,chemical_potentials!, rrate,ldebye
export showstruct

include("pnpsystem.jl")
export PNPSystem
export pnpunknowns,electrolytedata


include("cells.jl")
export voltagesweep,doublelayercap,bulkbc



end # module
