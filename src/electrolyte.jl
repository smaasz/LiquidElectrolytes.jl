

"""
$(TYPEDEF)

Data for electrolyte.

$(TYPEDFIELDS)
"""
@with_kw mutable struct ElectrolyteData{TBoundary}
    "Number of charged species."
    nc::Int=2

    "Potential index in species list."
    iϕ::Int=nc+1

    "Pressure index in species list"
    ip::Int=nc+2

    "Temperature"
    T::Float64=(273.15+25)*K
    
    "Mobility coefficient"
    D::Vector{Float64}=fill(2.0e-9m^2/s,nc)

    "Charge numbers"
    z::Vector{Int}=[ (-1)^(i-1) for i=1:nc]
    
    "Charge numbers scaled by F/RT"
    Z::Vector{Float64}=z.*F/(R*T)

    "Bulk concentration"
    c_bulk::Vector{Float64}=fill(0.1*mol/dm^3,nc)

    "Molar volume of solvent"
    v0::Float64=1/(55.4*mol/dm^3)

    "Solvation numbers"
    κ::Vector{Float64}=fill(10,nc)
    
    "Molar volumes"
    v::Vector{Float64}=fill(v0,nc)
    
    "Dielectric permittivity of water"
    ε::Float64=78.49

    "Boundary data"
    bdata::Union{Nothing,TBoundary}=nothing
end

ElectrolyteData(;kwargs...)=ElectrolyteData{Nothing}(;kwargs...)
ElectrolyteData(bdata;kwargs...)=ElectrolyteData{typeof(bdata)}(;bdata,kwargs...)
