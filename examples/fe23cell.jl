module fe23cell
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors, Parameters
using StaticArrays
using CompositeStructs
using StaticArrays

@phconstants N_A e R ε_0
const F=N_A*e

@siunits nm cm μF mol dm s mA A


@composite @kwdef mutable struct FE23Cell <: AbstractElectrolyteData
    ElectrolyteData...
    Γ_we::Int=1
    ϕ_we::Float64=0.0
    R0::Float64=1.0e-14*mol/(cm^2*s)
    Δg::Float64=0.0
    β::Float64=0.5
    ihplus::Int=1
    ife2::Int=2
    ife3::Int=3
    iso4::Int=4
end




function halfcellbc(f,u,bnode,data)
    bulkbc(f,u,bnode,data)
    if bnode.region==data.Γ_we
        @unpack R0,β,Δg,ihplus,ife2,ife3,iso4,iϕ,nc=data
        if !data.neutralflag
            boundary_dirichlet!(f,u,bnode;species=data.iϕ,region=data.Γ_we,value=data.ϕ_we)
        end
        μ0,μ=chemical_potentials!(MVector{4,eltype(u)}(undef),u,data)
        A=(μ[ife2]-μ[ife3]+Δg - data.neutralflag*F*(u[iϕ]-data.ϕ_we))/(R*data.T)
        r=rrate(R0,β,A)
        f[ife2]-=r
        f[ife3]+=r
    end
end


function main(;nref=0,compare=false,neutral=false,n=100,ϕmax=0.2, dlcap=false,R0=1.0e-6,logreg=1.0e-18,scheme=:μex,kwargs...)
    defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose=false,
              tol_relative=1.0e-14,
              tol_mono=1.0e-10)
    kwargs=merge(defaults, kwargs) 
    
    hmin=1.0e-1*nm*2.0^(-nref)
    hmax=1.0*nm*2.0^(-nref)
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)
    
    grid=simplexgrid(X)

    celldata=FE23Cell(;nc=4, z=[1,2,3,-2], neutralflag=neutral,κ=fill(0,4), Γ_we=1, Γ_bulk=2, R0=R0*mol/(cm^2*s),scheme,logreg)

    @unpack iϕ,ihplus,ife2,ife3,iso4, ip=celldata
    sc=1
    celldata.c_bulk[ihplus]=sc*1.0*mol/dm^3
    celldata.c_bulk[ife2]=sc*0.1*mol/dm^3
    celldata.c_bulk[ife3]=sc*0.1*mol/dm^3
    celldata.c_bulk[iso4]=sc*0.75*mol/dm^3

    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-12)
    @show celldata
    @show ldebye(celldata)
    
    cell=PNPSystem(grid;bcondition=halfcellbc,celldata)

    if compare
        celldata.neutralflag=false
        volts,currs, sols=voltagesweep(cell;ϕmax,n,ispec=ife2,kwargs...)
        
        tsol=VoronoiFVM.TransientSolution(sols,volts)
        
        for it=1:length(tsol.t)
            tsol.u[it][ife2,:]/=mol/dm^3
            tsol.u[it][ife3,:]/=mol/dm^3
        end
        
        celldata.neutralflag=true
        nvolts,ncurrs, sols=voltagesweep(cell;ϕmax,n,ispec=ife2,kwargs...)
        
        ntsol=VoronoiFVM.TransientSolution(sols,volts)
        
        vis=GridVisualizer(resolution=(600,400),Plotter=PyPlot,clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(A/m^2)")
        scalarplot!(vis,volts,-currs,color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP")
        scalarplot!(vis,nvolts,-ncurrs,clear=false,color=:green,markershape=:none,label="NNP")
        return reveal(vis)
    end
    
    #---------
    if dlcap
        molarities=[0.001,0.01,0.1,1]
    
    vis=GridVisualizer(resolution=(500,300),legend=:rt,clear=true,xlabel="φ/V",ylabel="C_dl/(μF/cm^2)",Plotter=PyPlot)
        hmol=1/length(molarities)
        for imol=1:length(molarities)
	    c=RGB(1-imol*hmol,0,imol*hmol)
	    t=@elapsed volts,caps=doublelayercap(cell;ϕmax,n,molarity=molarities[imol],kwargs...)
	    cdl0=Cdl0(cell.physics.data)
	    scalarplot!(vis,volts,caps/(μF/cm^2),
		        color=c,clear=false,label="$(molarities[imol])M",markershape=:none)
	    scalarplot!(vis,[0],[cdl0]/(μF/cm^2),clear=false,markershape=:circle,label="")
        end
        return
    end
    volts,currs, sols=voltagesweep(cell;ϕmax,n,ispec=ife2,kwargs...)
    tsol=VoronoiFVM.TransientSolution(sols,volts)

    for it=1:length(tsol.t)
        tsol.u[it][ife2,:]/=mol/dm^3
        tsol.u[it][ife3,:]/=mol/dm^3
    end

    @show extrema(tsol[ife2,end,:])
    @show extrema(tsol[ife3,end,:])
    if true
        xmax=0.25*nm
#        xmax=L
        xlimits=[0,xmax]
        vis=GridVisualizer(resolution=(1200,400),layout=(1,5),Plotter=PyPlot,clear=true)
        aspect=[2*xmax/(ϕmax)]
        scalarplot!(vis[1,1],F*currs/(mA/cm^2),volts,markershape=:none,title="IV",xlabel="I",ylabel="V")
        scalarplot!(vis[1,2],cell,tsol;species=ife2,aspect,xlimits,title="Fe2+",colormap=:summer)
        scalarplot!(vis[1,3],cell,tsol;species=ife3,aspect,xlimits,title="Fe3+",colormap=:summer)
        scalarplot!(vis[1,4],cell,tsol;species=iϕ,aspect,xlimits,title="ϕ",colormap=:bwr)
        scalarplot!(vis[1,5],cell,tsol;species=ip,aspect,xlimits,title="p",colormap=:summer)
    else
    
        vis=GridVisualizer(resolution=(500,300),Plotter=PyPlot,legend=:rb)
        scalarplot!(vis,grid, tsol[ife2,:,end], color=:red, label="fe2+")
#        scalarplot!(vis,grid, tsol[iϕ,:,1], color=:red, label="ϕ",flimits=(-2.5e-5,2.5e-5))
        scalarplot!(vis,grid, tsol[ife3,:,end], color=:green, label="fe3+",clear=false)
    
    end
#    save(plotsdir("1DResults.pdf"),vis)
    reveal(vis)
end

end
