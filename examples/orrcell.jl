module orrcell
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors, Parameters
using StaticArrays
using CompositeStructs

@phconstants R 
@siunits nm cm μF mol dm s


@composite @kwdef mutable struct ORRCell
    ElectrolyteData...
    TwoElectrodeCell...
end

function main(;n=100,ϕmax=0.2,molarity=1.0, kwargs...)
    defaults=(; max_round=3,tol_round=1.0e-9, verbose=true, tol_relative=1.0e-7,tol_mono=1.0e-10)
    kwargs=merge(defaults, kwargs) 
    
    ihplus::Int=1
    iso4::Int=2
    io2::Int=3
    
    
    celldata=ORRCell(;nc=3, z=[1,-2,0], κ=[0,0,0])
    celldata.c_bulk[io2]=0.001*mol/dm^3
    celldata.c_bulk[iso4]=molarity*mol/dm^3
    celldata.c_bulk[ihplus]=2.0*molarity*mol/dm^3




    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-12)

    
    R0::Float64=10.0e-13*mol/(cm^2*s)
    Δg::Float64=0.0
    β::Float64=0.5
    iϕ=celldata.iϕ
    
    function halfcellbc(f,u,bnode,data)
        bulkbc(f,u,bnode,data)
        if bnode.region==data.Γ_we
            μh2o,μ=chemical_potentials!(f,u,data)
            @show value(μh2o)
            A=-(4*μ[ihplus]+μ[io2]-2μh2o+Δg)/(R*data.T)
            r=rrate(R0,β,A)
            f.=0.0
            f[ihplus]+=4*r
            f[io2]+=r
            boundary_dirichlet!(f,u,bnode;species=data.iϕ,region=data.Γ_we,value=data.ϕ_we)
        end
    end
    
    hmin=1.0e-3*nm
    hmax=1.0*nm
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)


    grid=simplexgrid(X)
    cell=NPPSystem(grid;bcondition=halfcellbc,celldata)

    vis=GridVisualizer(resolution=(1200,400),layout=(1,5),Plotter=PyPlot)

    volts,currs, sols=voltagesweep(cell;ϕmax,n,ispec=io2,kwargs...)
    tsol=VoronoiFVM.TransientSolution(sols,volts)

    for it=1:length(tsol.t)
        tsol.u[it][io2,:]/=mol/dm^3
        tsol.u[it][ihplus,:]/=mol/dm^3
        tsol.u[it][iso4,:]/=mol/dm^3
    end

    
    xlimits=[0,2*nm]
    aspect=[2.5*nm/ϕmax]
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
