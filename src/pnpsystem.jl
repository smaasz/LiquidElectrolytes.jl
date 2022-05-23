
function pnpstorage(f, u, node, electrolyte)
    f[electrolyte.iϕ] = zero(eltype(u))
    f[electrolyte.ip] = zero(eltype(u))
    for ic = 1:electrolyte.nc
        f[ic] = u[ic]
    end
end


charge(u,i,electrolyte) = @views charge(u[:,i], electrolyte)


function pnpreaction(f, u, node, electrolyte)
    ## Charge density
    f[electrolyte.iϕ] = -charge(u,electrolyte)
    for ic = 1:electrolyte.nc
        f[ic] = 0
    end
end


default_bcondition(f,u,bnode,electrolyte)= nothing

function bulkbc(f,u,bnode,data)
    @unpack iϕ,ip,nc,Γ_bulk,ϕ_bulk,p_bulk,c_bulk=data
    if bnode.region==Γ_bulk
        boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_bulk,value=ϕ_bulk)
        boundary_dirichlet!(f,u,bnode,species=ip,region=Γ_bulk,value=p_bulk)
        for ic=1:nc
            boundary_dirichlet!(f,u,bnode,species=ic,region=Γ_bulk,value=data.c_bulk[ic])
        end
    end
end



"""
 Sedan flux

 Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

 see also the 198? Fortran code available via http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function pnpflux(f, u, edge, electrolyte)
    iϕ = electrolyte.iϕ # index of potential
    ip = electrolyte.ip
    pk = u[ip,1]
    pl = u[ip,2]
    ## Poisson flux
    dϕ = u[iϕ,1] - u[iϕ,2]
    f[iϕ]=electrolyte.ε*ε_0*dϕ*!electrolyte.neutralflag
    dp = pk-pl
    q1=charge(u,1,electrolyte)
    q2=charge(u,2,electrolyte)
    
    f[ip]=dp + (q1+q2)*dϕ/(2*electrolyte.pscale)

    if !iszero(electrolyte.v)
        # (log_c0k, log_bar_ck) = log_c0_barc(u, 1, electrolyte)
        # (log_c0l, log_bar_cl) = log_c0_barc(u, 2, electrolyte)


        c0k,bar_ck=c0_barc(u, 1, electrolyte)
        c0l,bar_cl=c0_barc(u, 2, electrolyte)
    end

    for ic = 1:electrolyte.nc
        if !iszero(electrolyte.v)
            M=electrolyte.M[ic]/electrolyte.M0
            V=electrolyte.v[ic]+(electrolyte.κ[ic]-M)*electrolyte.v0
            #            muk = V*pk/(R*electrolyte.T) -M*log_c0k + (M-1.0)*log_bar_ck
            #            mul = V*pl/(R*electrolyte.T) -M*log_c0l + (M-1.0)*log_bar_cl
            
            muk = V*pk +log(bar_ck^(M-1.0)/(c0k^M))*(R*electrolyte.T)
            mul = V*pl +log(bar_cl^(M-1.0)/(c0l^M))*(R*electrolyte.T)
        else
            muk=0.0
            mul=0.0
        end
        
        ## Combine potential gradient and excess chemical gradient
        arg = electrolyte.Z[ic] * (u[iϕ,1] - u[iϕ,2])  + (muk - mul)   /(R*electrolyte.T)

        ## Call Bernoulli function
        bp, bm = fbernoulli_pm(arg)

        ## Calculate drift-diffusion flux
        f[ic] = electrolyte.D[ic] * (bm * u[ic,1] - bp * u[ic,2])
    end
end 


"""
Averaging of reciprocal activity coefficients

"""
function pnpaflux(f, u, edge, electrolyte)
    iϕ = electrolyte.iϕ # index of potential
    ip = electrolyte.ip
    pk = u[ip,1]
    pl = u[ip,2]
    ## Poisson flux
    dϕ = u[iϕ,1] - u[iϕ,2]
    f[iϕ]=electrolyte.ε*ε_0*dϕ*!electrolyte.neutralflag
    dp = pk-pl
    q1=charge(u,1,electrolyte)
    q2=charge(u,2,electrolyte)
    
    f[ip]=dp + (q1+q2)*dϕ/(2*electrolyte.pscale)

    if !iszero(electrolyte.v)
        c0k,bar_ck=c0_barc(u, 1, electrolyte)
        c0l,bar_cl=c0_barc(u, 2, electrolyte)
    end
    
    for ic = 1:electrolyte.nc
        if !iszero(electrolyte.v)
            M=electrolyte.M[ic]/electrolyte.M0
            V=electrolyte.v[ic]+(electrolyte.κ[ic]-M)*electrolyte.v0
            γk=exp(V*pk/(R*electrolyte.T))*bar_ck^(M-1.0)/(c0k^M) 
            γl=exp(V*pl/(R*electrolyte.T))*bar_cl^(M-1.0)/(c0l^M) 
        else
            γk=1.0
            γl=1.0
        end

        ## Combine potential gradient and excess chemical gradient
        arg = electrolyte.Z[ic] * (u[iϕ,1] - u[iϕ,2])

        ## Call Bernoulli function
        bp, bm = fbernoulli_pm(arg)

        ## Calculate drift-diffusion flux
        f[ic] = electrolyte.D[ic] * (bm * u[ic,1]*γk - bp * u[ic,2]*γl)*2.0/(γk+γl)
    end
end 

function pnpcflux(f, u, edge, electrolyte)
    iϕ = electrolyte.iϕ # index of potential
    ip = electrolyte.ip
    pk = u[ip,1]
    pl = u[ip,2]
    ## Poisson flux
    dϕ = u[iϕ,1] - u[iϕ,2]
    f[iϕ]=electrolyte.ε*ε_0*dϕ*!electrolyte.neutralflag
    dp = pk-pl
    q1=charge(u,1,electrolyte)
    q2=charge(u,2,electrolyte)
    
    f[ip]=dp + (q1+q2)*dϕ/(2*electrolyte.pscale)

    if !iszero(electrolyte.v)
        c0k,bar_ck=c0_barc(u, 1, electrolyte)
        c0l,bar_cl=c0_barc(u, 2, electrolyte)
    end
    
    for ic = 1:electrolyte.nc
        ck = u[ic,1]
        cl = u[ic,2]
        if !iszero(electrolyte.v)
            M=electrolyte.M[ic]/electrolyte.M0
            V=electrolyte.v[ic]+(electrolyte.κ[ic]-M)*electrolyte.v0
            muk = V*pk +log(bar_ck^(M-1.0)/(c0k^M))*(R*electrolyte.T)
            mul = V*pl +log(bar_cl^(M-1.0)/(c0l^M))*(R*electrolyte.T)
        else
            muk=0.0
            mul=0.0
        end
        xlog(u)= u < 1.0e-20 ? log(1.0e-20) : log(u)

        hk = xlog(ck/bar_ck)*(R*electrolyte.T) - muk
        hl = xlog(cl/bar_cl)*(R*electrolyte.T) - mul

        f[ic] = electrolyte.D[ic] * 0.5 * (ck + cl) * (hk - hl + electrolyte.Z[ic] *(u[iϕ,1] - u[iϕ,2]))/(R*electrolyte.T)
    end
    @show value(f[1])
end 



function PNPSystem(grid;celldata=nothing,bcondition=default_bcondition,kwargs...)
    sys=VoronoiFVM.System(grid;
                          data=celldata,
                          flux=pnpflux,
                          reaction=pnpreaction,
                          storage=pnpstorage,
                          bcondition,
                          species=[ 1:celldata.nc..., celldata.iϕ,celldata.ip],
                          kwargs...
                          )
end

electrolytedata(sys)=sys.physics.data

function pnpunknowns(sys)
    @unpack iϕ,ip,nc,c_bulk=electrolytedata(sys)
    u=unknowns(sys)
    @views u[iϕ,:] .= 0
    @views u[ip,:] .= 0
    for ic=1:nc
        @views u[ic,:] .= c_bulk[ic]
    end
    u
end

