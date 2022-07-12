
"""
$(TYPEDEF)

Abstract super type for electrolytes
"""
abstract type AbstractElectrolyteData end


"""
$(TYPEDEF)

Data for electrolyte.

$(TYPEDFIELDS)
"""
@with_kw mutable struct ElectrolyteData <: AbstractElectrolyteData
    "Number of ionic species."
    nc::Int=2

    "Potential index in species list."
    iϕ::Int=nc+1

    "Pressure index in species list"
    ip::Int=nc+2

    "Temperature"
    T::Float64=(273.15+25)*K
    
    "Mobility coefficient"
    D::Vector{Float64}=fill(2.0e-9m^2/s,nc)
    
    "Molar weight of solvent"
    M0::Float64=18.0153*g/mol

    "Molar weight of ions"
    M::Vector{Float64}=fill(M0,nc)

    "Charge numbers of ions"
    z::Vector{Int}=[ (-1)^(i-1) for i=1:nc]
    
    "Charge numbers scaled by F/RT"
    Z::Vector{Float64}=z.*F/(R*T)

    "Bulk ion concentrations"
    c_bulk::Vector{Float64}=fill(0.1*mol/dm^3,nc)

    "Bulk voltage"
    ϕ_bulk::Float64=0.0
    
    "Bulk pressure"
    p_bulk::Float64=0.0

    "Bulk boundary number"
    Γ_bulk::Int=1
    
    "Molar volume of solvent"
    v0::Float64=1/(55.4*mol/dm^3)

   "Solvation numbers"
    κ::Vector{Float64}=fill(0,nc)
    
    "Molar volumes of ions"
    v::Vector{Float64}=fill(v0,nc)
    
    "Dielectric permittivity of solvent"
    ε::Float64=78.49

    "Pressure scaling factor"
    pscale::Float64=1.0e9

    "Electroneutrality switch"
    neutralflag::Bool=false

    "Electroneutrality switch"
    scheme::Symbol=:μex

    "Regularization parameter"
    epsreg::Float64=1.0e-20
end

Base.show(io::IO, this::AbstractElectrolyteData)=showstruct(io,this)

"""
    Cdl0(electrolyte)

Double layer capacitance at zero voltage for symmetric binary electrolyte.

### Example

```jldoctest
using LessUnitful
@unitfactors mol dm μF cm
ely=ElectrolyteData(c_bulk=fill(0.01*mol/dm^3,2))
round(Cdl0(ely)/(μF/cm^2),digits=2)
# output

22.85
```
"""
Cdl0(data::AbstractElectrolyteData)=sqrt( 2*(data.ε)*ε_0*F^2*data.c_bulk[1]/(R*data.T));


"""
    ldebye(electrolyte)

Debye length.
"""
ldebye(data)=sqrt( data.ε*ε_0*R*data.T/(F^2*data.c_bulk[1]))



"""
    charge(c,electrolyte)

Calculate charge from vector of concentrations
"""
function charge(u,electrolyte::AbstractElectrolyteData)
    q=zero(eltype(u))
    for ic=1:electrolyte.nc
        q+=u[ic] * electrolyte.z[ic]
    end
    q*F
end

@doc raw"""
	vrel(ic,electrolyte)

``v_{i,rel}=κ_i+\frac{v_i}{v_0}``
"""
vrel(ic,electrolyte)=electrolyte.v[ic]/electrolyte.v0+electrolyte.κ[ic]


@doc raw"""
	c0_barc(u,electrolyte)

Calculate solvent concentration ``c_0`` and summary concentration ``\bar c`` from vector of concentrations `c`
using the incompressibility constraint (assuming ``κ_0=0``):
```math
 \sum_{i=0}^N c_i (v_i + κ_iv_0) =1
```

This gives

```math
 c_0v_0=1-\sum_{i=1}^N c_i (v_i+ κ_iv_0)
```

```math
c_0= 1/v_0 - \sum_{i=1}^N c_i(\frac{v_i}{v_0}+κ)
```

Then we can calculate 
```math
 \bar c= \sum_{i=0}^N c_i
```
"""
function c0_barc(c, electrolyte)
    c0 = one(eltype(c)) / electrolyte.v0
    barc = zero(eltype(c))
    for ic = 1:electrolyte.nc
        barc += c[ic]
        c0 -= c[ic] * vrel(ic,electrolyte)
    end
    barc += c0
    c0+electrolyte.epsreg, barc+electrolyte.epsreg
end

"""
    rlog(u, electrolyte)

Regularized logarithm:

```
   rlog(u,electrolyte)= log(u+electrolyte.epsreg)
```
"""
rlog(x,electrolyte::AbstractElectrolyteData)=rlog(x,eps=electrolyte.epsreg)

function rlog(x;eps=1.0e-20)
    if x<eps
        return log(eps)+(x-eps)/eps
    else
        return log(x)
    end
end



function c0(U::Array, electrolyte)
    c0 = similar(U[1,:])
    c0 .= 1.0 / electrolyte.v0
    for ic = 1:electrolyte.nc
        c0 -= U[ic,:] .* vrel(ic,electrolyte)
    end
    c0
end

"""
    chemical_potentials!(μ,u,electrolyte::AbstractElectrolyteData)

Calculate chemical potentials from concentrations.

Input:
  -  `μ`: mutated, allocated memory for result
  -  `u`: concentrations
Returns `μ0, μ`: chemical potential of solvent and chemical
potentials of ions.

"""
function chemical_potentials!(μ,u,data::AbstractElectrolyteData)
    c0,barc=c0_barc(u,data)
    p=u[data.ip]*data.pscale
    p_ref=0
    μ0=rlog(c0/barc,data)*R*data.T+data.v0*(p-p_ref)
    for i=1:data.nc
        μ[i]=rlog(u[i]/barc,data)*R*data.T+data.v[i]*(p-p_ref)
    end
    μ0,μ
end


"""
    rrate(R0,β,A)

Thermodynamic reaction rate expression


    rrate(R0,β,A)=R0*(exp(-β*A) - exp((1-β)*A))
"""
rrate(R0,β,A)=R0*(exp(-β*A) - exp((1-β)*A))

"""
    wnorm(u,w,p)

Weighted norm with respect to columns
"""
function wnorm(u,w,p)
    @views norms=[w[i]*LinearAlgebra.norm(u[i,:],p) for i=1:size(u,1)]
    LinearAlgebra.norm(norms,p)
end
