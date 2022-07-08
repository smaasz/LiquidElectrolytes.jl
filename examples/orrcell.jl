module orrcell
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors, Parameters
using StaticArrays
using LessUnitful
using CompositeStructs

@phconstants R N_A e
@siunits nm cm μF mol dm s
const F=N_A*e


@composite @kwdef mutable struct ORRCell <: AbstractElectrolyteData
    ElectrolyteData...
    R0::Float64=10.0e-6*mol/(cm^2*s)
    Δg::Float64=0.0
    β::Float64=0.5
    Γ_we::Int=1
    ϕ_we::Float64=0.0
    ihplus::Int=1
    iso4::Int=2
    io2::Int=3
end



function halfcellbc(f,u,bnode,data)
    bulkbc(f,u,bnode,data)
    if bnode.region==data.Γ_we
        @unpack R0,β,Δg,iϕ,ihplus,iso4,io2=data
        f.=0.0
        if !data.neutralflag
            boundary_dirichlet!(f,u,bnode;species=data.iϕ,region=data.Γ_we,value=data.ϕ_we)
        end
        μh2o,μ=chemical_potentials!(MVector{4,eltype(u)}(undef),u,data)
        A=(4*μ[ihplus]+μ[io2]-2μh2o+Δg + data.neutralflag*F*(u[iϕ]-data.ϕ_we))/(R*data.T)
        r=rrate(R0,β,A)
        f[ihplus]-=4*r
        f[io2]-=r
    end
end


function main(;voltages=-0.2:0.005:0.2,molarity=0.1,nref=0,neutral=false,scheme=:μex,logreg=1.0e-10,kwargs...)
    defaults=(; max_round=3,tol_round=1.0e-10, verbose=true, tol_relative=1.0e-7,tol_mono=1.0e-10)
    kwargs=merge(defaults, kwargs) 
    
    hmin=1.0e-1*nm*2.0^(-nref)
    hmax=1.0*nm*2.0^(-nref)
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)

    grid=simplexgrid(X)
    
    
    celldata=ORRCell(;nc=3, z=[1,-2,0], κ=fill(0,3), Γ_we=1, Γ_bulk=2,neutralflag=neutral,logreg,scheme)

    @unpack iϕ,ihplus,iso4,io2=celldata

    
    celldata.c_bulk[io2]=0.001*mol/dm^3
    celldata.c_bulk[iso4]=molarity*mol/dm^3
    celldata.c_bulk[ihplus]=2.0*molarity*mol/dm^3
    celldata.R0=10.0e-8*mol/(cm^2*s)


    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-12)
    
    cell=PNPSystem(grid;bcondition=halfcellbc,celldata)
    check_allocs!(cell,false)
    
    vis=GridVisualizer(resolution=(1200,400),layout=(1,5),Plotter=PyPlot)

    volts,currs, sols=voltagesweep(cell;voltages,ispec=io2,kwargs...)
    tsol=VoronoiFVM.TransientSolution(sols,volts)

    for it=1:length(tsol.t)
        tsol.u[it][io2,:]/=mol/dm^3
        tsol.u[it][ihplus,:]/=mol/dm^3
        tsol.u[it][iso4,:]/=mol/dm^3
    end

    xmax=20*nm
    xlimits=[0,xmax]
        aspect=3.5*xmax/(tsol.t[end]-tsol.t[begin])

    scalarplot!(vis[1,1],currs,volts,markershape=:none,title="IV",xlabel="I",ylabel="V")
    scalarplot!(vis[1,2],cell,tsol;species=io2,aspect,xlimits,title="O2",colormap=:summer)
    scalarplot!(vis[1,3],cell,tsol;species=ihplus,aspect,xlimits,title="H+",colormap=:summer)
    scalarplot!(vis[1,4],cell,tsol;species=iϕ,aspect,xlimits,title="ϕ",colormap=:bwr)
    scalarplot!(vis[1,5],tsol[io2,1,:],volts,title="c_o2(0)",xlabel="O2",ylabel="V")

#    Plots.plot(tsol.t,tsol[io2,1,:])
   
#    save(plotsdir("1DResults.pdf"),vis)
    reveal(vis)
end

end
