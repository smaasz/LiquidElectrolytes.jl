"""
$(TYPEDEF)

Abstract super type for electrolytes
"""
abstract type AbstractElectrolyteData end


"""
$(TYPEDEF)

Data for electrolyte. It is defined using [Parameters.jl](https://github.com/mauro3/Parameters.jl), therefore it
has keyword constructors like
```julia
    ElectrolyteData(nc=3,z=[-1,2,1])
```

$(TYPEDFIELDS)
"""
@kwdef mutable struct ElectrolyteData <: AbstractElectrolyteData
    "Number of ionic species."
    nc::Int=2
    
    "Potential index in species list."
    iϕ::Int=nc+1
    
    "Pressure index in species list"
    ip::Int=nc+2
    
    "Mobility coefficient"
    D::Vector{Float64}=fill(2.0e-9*ufac"m^2/s",nc) 

    "Charge numbers of ions"
    z::Vector{Int}=[ (-1)^(i-1) for i=1:nc]
    
    "Molar weight of solvent"
    M0::Float64=18.0153*ufac"g/mol"
    
    "Molar weight of ions"
    M::Vector{Float64}=fill(M0,nc)

    "Molar volume of solvent"
    v0::Float64=1/(55.4*ufac"M")

    "Molar volumes of ions"
    v::Vector{Float64}=fill(v0,nc)
    
    "Solvation numbers"
    κ::Vector{Float64}=fill(10.0,nc)
    
    "Bulk ion concentrations"
    c_bulk::Vector{Float64}=fill(0.1*ufac"M",nc)
    
    "Bulk voltage"
    ϕ_bulk::Float64=0.0*ufac"V"
    
    "Bulk pressure"
    p_bulk::Float64=0.0*ufac"Pa"

    "Bulk boundary number"
    Γ_bulk::Int=2

    "Working electrode voltage"
    ϕ_we::Float64=0.0*ufac"V"
    
    "Working electrode  boundary number"
    Γ_we::Int=1
    
    "Temperature"
    T::Float64=(273.15+25)*ufac"K"
    
    "Molar gas constant scaled with temperature"
    RT::Float64=ph"R"*T
    
    "Faraday constant"
    F::Float64=ph"N_A*e"
    
    "Dielectric permittivity of solvent"
    ε::Float64=78.49
    
    "Dielectric permittivity of vacuum"
    ε_0::Float64=ph"ε_0"
    
    "Pressure scaling factor"
    pscale::Float64=1.0e9
    
    "Local electroneutrality switch"
    eneutral::Bool=false
    
    """
    [Flux caculation scheme](@id fluxes)
    This allows to choose between
    - `:μex` (default): excess chemical potential (SEDAN) scheme, see [`sflux`](@ref)
    - `:act` : scheme based on reciprocal activity coefficients, see [`aflux`](@ref)
    - `:cent` : central scheme, see [`cflux`](@ref).
    """
    scheme::Symbol=:μex
    
    """
    Regularization parameter used in [`rlog`](@ref)
    """
    epsreg::Float64=1.0e-20
end

function Base.show(io::IO, ::MIME"text/plain", this::ElectrolyteData)
    showstruct(io,this)
end

"""
    dlcap0(electrolyte)

Double layer capacitance at zero voltage for symmetric binary electrolyte.

### Example
```jldoctest
using LessUnitful
ely = ElectrolyteData(c_bulk=fill(0.01ufac"mol/dm^3",2))
round(dlcap0(ely),sigdigits=5) |> u"μF/cm^2"
# output

22.847 μF cm^-2
```
"""
function dlcap0(data::AbstractElectrolyteData)
    sqrt( 2*(data.ε)*data.ε_0*data.F^2*data.c_bulk[1]/(data.RT))
end


"""
    debyelength(electrolyte)

Debye length.

```jldoctest
using LessUnitful
ely = ElectrolyteData(c_bulk=fill(0.01ufac"mol/dm^3",2))
round(debyelength(ely),sigdigits=5) |> u"nm"
# output

4.3018 nm
```
"""
debyelength(data)=sqrt(data.ε*data.ε_0*data.RT/(data.F^2*data.c_bulk[1]))



"""
    charge(c,electrolyte)

Calculate charge from vector of concentrations
"""
function charge(u,electrolyte::AbstractElectrolyteData)
    q=zero(eltype(u))
    for ic=1:electrolyte.nc
        q+=u[ic] * electrolyte.z[ic]
    end
    q*electrolyte.F
end

@doc raw"""
	vrel(ic,electrolyte)

Calculate relative (wrt. solvent) molar volume of i-th species ``v_{i,rel}=κ_i+\frac{v_i}{v_0}``.
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

Calls rlog(u;eps=electrolyte.epsreg)
"""
rlog(x,electrolyte::AbstractElectrolyteData)=rlog(x,eps=electrolyte.epsreg)

"""
    rlog(u; eps=1.0e-20)

Regularized logarithm. Smooth linear continuation for `x<eps`.
This means we can calculate a "logarithm"  of a small negative number.
"""
function rlog(x;eps=1.0e-20)
    if x<eps
        return log(eps)+(x-eps)/eps
    else
        return log(x)
    end
end

"""
       solventconcentration(U::Array, electrolyte)

Calculate vector of solvent concentrations from solution array.
"""
function solventconcentration(U::Array, electrolyte)
    c0 = similar(U[1,:])
    c0 .= 1.0 / electrolyte.v0 + electrolyte.epsreg
    for ic = 1:electrolyte.nc
        c0 -= U[ic,:] .* vrel(ic,electrolyte)
    end
    c0
end

"""
    chemical_potentials!(μ,u,electrolyte)

Calculate chemical potentials from concentrations.

Input:
  -  `μ`: memory for result (will be filled)
  -  `u`: local solution vector (concentrations, potential, pressure)
Returns `μ0, μ`: chemical potential of solvent and chemical potentials of ions.


```jldoctest
using LessUnitful
ely = ElectrolyteData(c_bulk=fill(0.01ufac"mol/dm^3",2))
μ0,μ=chemical_potentials!([0.0,0.0],vcat(ely.c_bulk,[0,0]),ely)
round(μ0,sigdigits=5),round.(μ,sigdigits=5)
# output

(-0.89834, [-21359.0, -21359.0])
```
"""
function chemical_potentials!(μ,u,data::AbstractElectrolyteData)
    c0,barc=c0_barc(u,data)
    p=u[data.ip]*data.pscale
    p_ref=0
    μ0=rlog(c0/barc,data)*data.RT+data.v0*(p-p_ref)
    for i=1:data.nc
        μ[i]=rlog(u[i]/barc,data)*data.RT+data.v[i]*(p-p_ref)
    end
    μ0,μ
end


"""
    rexp(x;trunc=500.0)

Regularized exponential. Linear continuation for `x>trunc`,  
returns 1/rexp(-x) for `x<trunc`.
"""
function rexp(x;trunc=500.0)
    if x<-trunc
        1.0/rexp(-x;trunc)
    elseif x<=trunc
        exp(x)
    else
        exp(trunc)*(x-trunc+1)
    end
end


"""
    rrate(R0,β,A)

Reaction rate expression

    rrate(R0,β,A)=R0*(exp(-β*A) - exp((1-β)*A))
"""
rrate(R0,β,A)=R0*(rexp(-β*A) - rexp((1-β)*A))

"""
    wnorm(u,w,p)

Weighted norm with respect to columns
"""
function wnorm(u,w,p)
    @views norms=[w[i]*LinearAlgebra.norm(u[i,:],p) for i=1:size(u,1)]
    LinearAlgebra.norm(norms,p)
end
