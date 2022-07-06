
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
    f[electrolyte.ip] = 0
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

function dμ(γk, γl, electrolyte)
    if abs(γk-γl)<1.0e-20
        return 0.0
    elseif abs(γk)<abs(γl)
        return -rlog(γk/γl,electrolyte)*(R*electrolyte.T)
    elseif abs(γk)>abs(γl)
        return rlog(γl/γk,electrolyte)*(R*electrolyte.T)
    else
        return 0.0
    end
end

function sflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)
    bp, bm = fbernoulli_pm(electrolyte.Z[ic] * dϕ  + dμ(γk,γl,electrolyte) /(R*electrolyte.T))
    electrolyte.D[ic] * (bm*ck - bp*cl)
end

function aflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)
 #   γk+=electrolyte.logreg
 #   γl+=electrolyte.logreg
    bp, bm = fbernoulli_pm(electrolyte.Z[ic] * dϕ)
    electrolyte.D[ic] * (bm*ck/γk - bp*cl/γl)*(γk+γl)/2
end

function cflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte)
    lck = rlog(ck/bar_ck,electrolyte)*(R*electrolyte.T)
    lcl = rlog(cl/bar_cl,electrolyte)*(R*electrolyte.T)
    electrolyte.D[ic] * 0.5 * (ck + cl) * (lck - lcl +  dμ(γk,γl,electrolyte)  + electrolyte.z[ic]*F*dϕ)/(R*electrolyte.T)
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
    @unpack ip, iϕ, v0, v, M0, M, κ,  ε, T, nc, neutralflag, pscale = electrolyte
    
    pk,pl = u[ip,1],u[ip,2]
    ϕk,ϕl = u[iϕ,1],u[iϕ,2]

    qk,ql=charge(u,1,electrolyte),charge(u,2,electrolyte)
    @views c0k,bar_ck=c0_barc(u[:,1], electrolyte)
    @views c0l,bar_cl=c0_barc(u[:,2], electrolyte)
    
    dϕ = ϕk-ϕl
    dp = pk-pl


    f[iϕ]=ε*ε_0*dϕ*!neutralflag
    f[ip]=dp + (qk+ql)*dϕ/(2*pscale)

    γk,γl=1.0,1.0
    

    for ic = 1:nc
        f[ic]=0.0
        ck,cl=u[ic,1],u[ic,2]
        if !iszero(v)
            Mrel=M[ic]/M0
            V=v[ic]+(κ[ic]-Mrel)*v0
            γk = exp(-V*pk/(R*T))*c0k^Mrel/bar_ck^(Mrel-1.0)
            γl = exp(-V*pl/(R*T))*c0l^Mrel/bar_cl^(Mrel-1.0)
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





function PNPSystem(grid;celldata=ElectrolyteData(),bcondition=default_bcondition,kwargs...)
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

