
function nppstorage(f, u, node, data)
    f[data.iϕ] = zero(eltype(u))
    for ic = 1:data.nc
        f[ic] = u[ic]
    end
end

function charge(c, data)
    q=zero(eltype(c))
    for ic=1:data.nc
        q+=c[ic] * data.z[ic]
    end
    q*F
end

charge(u,i,data) = @views charge(u[:,i], data)


function nppreaction(f, u, node, data)
    ## Charge density
    f[data.iϕ] = -charge(u,data)
    for ic = 1:data.nc
        f[ic] = 0
    end
end


default_bcondition(f,u,bnode,data)= nothing


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

vrel(ic,data)=data.v[ic]/data.v0+data.κ[ic]
    
function c0_barc(c, data)
    c0 = one(eltype(c)) / data.v0
    barc = zero(eltype(c))
    for ic = 1:data.nc
        barc += c[ic]
        c0 -= c[ic] * vrel(ic,data)
    end
    barc += c0
    c0, barc
end

xlog(u)= u<1.0e-20 ? -20.0*one(u) : log(u)
#xlog(u)=  log(u)

c0_barc(u,i,data) = @views c0_barc(u[:,i], data)

log_c0_barc(u,i,data) = @views xlog.(c0_barc(u, i, data))

function csol(U::Array, data)
    C0 = similar(U[1,:])
    C0 .= 1.0 / data.v0
    for ic = 1:data.nc
        C0 -= U[ic,:] .* vrel(ic,data)
    end
    C0
end

"""
 Sedan flux

 Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

 see also the 198? Fortran code available via http://www-tcad.stanford.edu/tcad/programs/oldftpable.html
"""
function nppflux(f, u, edge, data)
    iϕ = data.iϕ # index of potential
    ip = data.ip
    ## Poisson flux
    dϕ = u[iϕ,1] - u[iϕ,2]
    dp = u[ip,1] - u[ip,2]

    q1=charge(u,1,data)
    q2=charge(u,2,data)

    f[iϕ]=data.ε*ε_0*dϕ

    f[ip]=dp + (q1+q2)*dϕ/2

    if !iszero(data.v)
        (log_c0k, log_bar_ck) = log_c0_barc(u, 1, data)
        (log_c0l, log_bar_cl) = log_c0_barc(u, 2, data)
    end
    
    for ic = 1:data.nc
        if !iszero(data.v)
            ## Caclculate excess chemical potentials
            muk = -log_bar_ck - (log_c0k - log_bar_ck) * vrel(ic,data)
            mul = -log_bar_cl - (log_c0l - log_bar_cl) * vrel(ic,data)
        else
            muk=0.0
            mul=0.0
        end
        
        ## Combine potential gradient and excess chemical gradient
        arg = data.Z[ic] * (u[iϕ,1] - u[iϕ,2]) + (muk - mul)

        ## Call Bernoulli function
        bp, bm = fbernoulli_pm(arg)
        
        ## Calculate drift-diffusion flux
        f[ic] = data.D[ic] * (bm * u[ic,1] - bp * u[ic,2])
    end
end 


function NPPSystem(grid;electrolyte=nothing,bcondition=default_bcondition,kwargs...)
    sys=VoronoiFVM.System(grid;
                          data=electrolyte,
                          flux=nppflux,
                          reaction=nppreaction,
                          storage=nppstorage,
                          bcondition,
                          species=[ 1:electrolyte.nc..., electrolyte.iϕ,electrolyte.ip],
                          kwargs...
                          )
end

electrolytedata(sys)::ElectrolyteData=sys.physics.data

boundarydata(sys)=sys.physics.data.bdata

function nppunknowns(sys)
    @unpack iϕ,ip,nc,c_bulk=electrolytedata(sys)
    u=unknowns(sys)
    @views u[iϕ,:] .= 0
    @views u[ip,:] .= 0
    for ic=1:nc
        @views u[ic,:] .= c_bulk[ic]
    end
    u
end

