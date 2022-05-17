
"""
$(TYPEDEF)

Data for electrolyte.

$(TYPEDFIELDS)
"""
@composite @Base.kwdef mutable struct ElectrolyteData
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
    
    "Molar weight of solvent"
    M0::Float64=18.0153*g/mol

    "Molar weight"
    M::Vector{Float64}=fill(M0,nc)

    "Charge numbers"
    z::Vector{Int}=[ (-1)^(i-1) for i=1:nc]
    
    "Charge numbers scaled by F/RT"
    Z::Vector{Float64}=z.*F/(R*T)

    "Bulk concentration"
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
    
    "Molar volumes"
    v::Vector{Float64}=fill(v0,nc)
    
    "Dielectric permittivity of water"
    ε::Float64=78.49


    "Pressure scaling factor"
    pscale::Float64=1.0e9
    
end


Cdl0(data)=sqrt( 2*(data.ε)*ε_0*F^2*data.c_bulk[1]/(R*data.T));

function charge(u,data)
    q=zero(eltype(u))
    for ic=1:data.nc
        q+=u[ic] * data.z[ic]
    end
    q*F
end

"""
Calculate c0 and \bar c
from using the incompressibility constraint
```math
 \\sum_{i=0}^N c_i v_i =1
```

This gives

```math
 c_0v_0=1-\\sum_{i=1}^N c_i v_i
 c_0= 1/v_0 - \\sum_{i=1}^N c_iv_i/v0
```

Then we can calculate 
```math
 \\bar c= \\sum_{i=0}^N c_i
```
"""

vrel(ic,electrolyte)=electrolyte.v[ic]/electrolyte.v0+electrolyte.κ[ic]
    
function c0_barc(c, electrolyte)
    c0 = one(eltype(c)) / electrolyte.v0
    barc = zero(eltype(c))
    for ic = 1:electrolyte.nc
        barc += c[ic]
        c0 -= c[ic] * vrel(ic,electrolyte)
    end
    barc += c0
    max(1.0e-20,c0), barc
end

xlog(u)= u<1.0e-50 ? -50.0*one(u) : log(u)
#xlog(u)=  log(u)

c0_barc(u,i,electrolyte) = @views c0_barc(u[:,i], electrolyte)

log_c0_barc(u,i,electrolyte) = @views log.(c0_barc(u, i, electrolyte))

function c0(U::Array, electrolyte)
    c0 = similar(U[1,:])
    c0 .= 1.0 / electrolyte.v0
    for ic = 1:electrolyte.nc
        c0 -= U[ic,:] .* vrel(ic,electrolyte)
    end
    max(c0,0.0)
end


function chemical_potentials!(μ,u,data)
    c0,barc=c0_barc(u,data)
    p=u[data.ip]*data.pscale
    p_ref=0
    μ0=xlog(c0/barc)*(R*data.T)+data.v0*(p-p_ref)
    for i=1:data.nc
        μ[i]=xlog(u[i]/barc)*(R*data.T)+data.v[i]*(p-p_ref)
    end
    μ0,μ
end

rrate(R0,β,A)=R0*(exp(-β*A) - exp((1-β)*A))

ldebye(data)=sqrt( data.ε*ε_0*R*data.T/(F^2*data.c_bulk[1]))

function wnorm(u,w,p)
    @views norms=[w[i]*LinearAlgebra.norm(u[i,:],p) for i=1:size(u,1)]
    LinearAlgebra.norm(norms,p)
end
