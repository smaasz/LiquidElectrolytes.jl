
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
        (log_c0k, log_bar_ck) = log_c0_barc(u, 1, electrolyte)
        (log_c0l, log_bar_cl) = log_c0_barc(u, 2, electrolyte)
    end
    
    for ic = 1:electrolyte.nc
        if !iszero(electrolyte.v)
            M=electrolyte.M[ic]/electrolyte.M0
            V=electrolyte.v[ic]+(electrolyte.κ[ic]-M)*electrolyte.v0
            muk = V*pk/(R*electrolyte.T) -M*log_c0k + (M-1.0)*log_bar_ck
            mul = V*pl/(R*electrolyte.T) -M*log_c0l + (M-1.0)*log_bar_cl
        else
            muk=0.0
            mul=0.0
        end

        ## Combine potential gradient and excess chemical gradient
        arg = electrolyte.Z[ic] * (u[iϕ,1] - u[iϕ,2]) + (muk - mul)

        ## Call Bernoulli function
        bp, bm = fbernoulli_pm(arg)

        ## Calculate drift-diffusion flux
        f[ic] = electrolyte.D[ic] * (bm * u[ic,1] - bp * u[ic,2])
    end
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

