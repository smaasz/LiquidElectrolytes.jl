
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

function sflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)
    μk=rlog(γk,electrolyte)*(R*electrolyte.T)
    μl=rlog(γl,electrolyte)*(R*electrolyte.T)
    bp, bm = fbernoulli_pm(electrolyte.Z[ic] * dϕ  + (μk - μl)   /(R*electrolyte.T))
    electrolyte.D[ic] * (bm*ck - bp*cl)
end

function aflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)
    bp, bm = fbernoulli_pm(electrolyte.Z[ic] * dϕ)
    electrolyte.D[ic] * (bm*ck*γk - bp*cl*γl)*2.0/(γk+γl)
end

function cflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)
    μk=rlog(γk,electrolyte)*(R*electrolyte.T)
    μl=rlog(γl,electrolyte)*(R*electrolyte.T)
    
    hk = rlog(ck/bar_ck,electrolyte)*(R*electrolyte.T) + μk
    hl = rlog(cl/bar_cl,electrolyte)*(R*electrolyte.T) + μl
    
    electrolyte.D[ic] * 0.5 * (ck + cl) * (hk - hl + electrolyte.z[ic]*F*dϕ)/(R*electrolyte.T)
end


"""
 Sedan flux

 Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

 see also the 198? Fortran code available via
  https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function pnpflux(f, u, edge, electrolyte)
    iϕ = electrolyte.iϕ # index of potential
    ip = electrolyte.ip

    pk,pl = u[ip,1],u[ip,2]
    ϕk,ϕl = u[iϕ,1],u[iϕ,2]

    qk,ql=charge(u,1,electrolyte),charge(u,2,electrolyte)

    c0k,bar_ck=c0_barc(u, 1, electrolyte)
    c0l,bar_cl=c0_barc(u, 2, electrolyte)
    
    dϕ = ϕk-ϕl
    dp = pk-pl

    f[iϕ]=electrolyte.ε*ε_0*dϕ*!electrolyte.neutralflag
    f[ip]=dp + (qk+ql)*dϕ/(2*electrolyte.pscale)

    γk,γl=1.0,1.0
    
    for ic = 1:electrolyte.nc
        ck,cl=u[ic,1],u[ic,2]
        if !iszero(electrolyte.v)
            M=electrolyte.M[ic]/electrolyte.M0
            V=electrolyte.v[ic]+(electrolyte.κ[ic]-M)*electrolyte.v0
            γk = exp(V*pk/(R*electrolyte.T))*bar_ck^(M-1.0)/(c0k^M)
            γl = exp(V*pl/(R*electrolyte.T))*bar_cl^(M-1.0)/(c0l^M)
        end
        
        if electrolyte.scheme==:μex
            f[ic]=sflux(ic,dϕ,ck,cl,γk,γl,bar_ck, bar_cl,electrolyte)
        elseif electrolyte.scheme==:act
            f[ic]=aflux(ic,dϕ,ck,cl,γk,γl,bar_ck, bar_cl,electrolyte)
        elseif  electrolyte.scheme==:cent
            f[ic]=cflux(ic,dϕ,ck,cl,γk,γl,bar_ck, bar_cl,electrolyte)
        else
            error("no such scheme: $(electrolyte.scheme)")
        end
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

