module LiquidElectrolytes

using VoronoiFVM,ExtendableGrids
using DocStringExtensions
using Parameters

include("units.jl")

@phconstants N_A e R ε_0
const F=N_A*e
@siunits K mol dm m s


export @siunits,@phconstants

include("electrolyte.jl")
include("halfcell.jl")
export ElectrolyteData, TBoundary
export ideally_polarizing_halfcell,oxygen_halfcell,voltagesweep

include("nppsystem.jl")
export NPPSystem
export nppunknowns,electrolytedata
export NPPHalfCell,doublelayercap

end # module
