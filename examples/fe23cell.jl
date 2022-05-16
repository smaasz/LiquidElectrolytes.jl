module fe23cell
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors, Parameters
using StaticArrays

@phconstants R 
@siunits nm cm μF mol dm s


function main(;n=100,ϕmax=0.15,kwargs...)
    defaults=(; max_round=3,tol_round=1.0e-9, verbose=false, tol_relative=1.0e-7,tol_mono=1.0e-10)
    kwargs=merge(defaults, kwargs) 
    
    ihplus::Int=1
    ife2::Int=2
    ife3::Int=3
    iso4::Int=4
    
    
    electrolyte=Electrolyte(;nc=4, z=[1,2,3,-2], κ=[0,0,0,0])
    electrolyte.c_bulk[ihplus]=1.0*mol/dm^3
    electrolyte.c_bulk[ife2]=0.1*mol/dm^3
    electrolyte.c_bulk[ife3]=0.1*mol/dm^3
    electrolyte.c_bulk[iso4]=0.75*mol/dm^3

    @assert isapprox(electrolyte.c_bulk'*electrolyte.z,0, atol=1.0e-12)
    

    @info electrolyte

    
    R0::Float64=1.0e-14*mol/(cm^2*s)
    Δg::Float64=0.0
    β::Float64=0.5
    iϕ=electrolyte.iϕ
    
    function halfcellbc(f,u,bnode,data)
        bd=data.bdata
        bulkbc(f,u,bnode,data)
        if bnode.region==bd.Γ_we
            μ0,μ=chemical_potentials!(f,u,data)
            A=(μ[ife2]-μ[ife3]+Δg)/(R*data.T)
            r=rrate(R0,β,A)
            f.=0.0
            f[ife2]+=r
            f[ife3]-=r
            boundary_dirichlet!(f,u,bnode;species=data.iϕ,region=bd.Γ_we,value=bd.ϕ_we)
        end
    end
    
    hmin=1.0e-3*nm
    hmax=1.0*nm
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)


    grid=simplexgrid(X)
    cell=NPPSystem(grid;bcondition=halfcellbc,electrolyte)


    volts,currs, sols=voltagesweep(cell;ϕmax,n,ispec=ife2,kwargs...)
    tsol=VoronoiFVM.TransientSolution(sols,volts)
    for it=1:length(tsol.t)
        tsol.u[it][ife2,:]/=mol/dm^3
        tsol.u[it][ife3,:]/=mol/dm^3
    end

    @show extrema(tsol[ife2,:,:])
    @show extrema(tsol[ife3,:,:])

    xlimits=[0,2*nm]
    vis=GridVisualizer(resolution=(1200,400),layout=(1,4),Plotter=PyPlot)
    aspect=[2.5*nm/ϕmax]
    scalarplot!(vis[1,1],currs,volts,markershape=:none,title="IV",xlabel="I",ylabel="V")
    scalarplot!(vis[1,2],cell,tsol;species=ife2,aspect,xlimits,title="Fe2+",colormap=:summer)
    scalarplot!(vis[1,3],cell,tsol;species=ife3,aspect,xlimits,title="Fe3+",colormap=:summer)
    scalarplot!(vis[1,4],cell,tsol;species=iϕ,aspect,xlimits,title="ϕ",colormap=:bwr)

#    save(plotsdir("1DResults.pdf"),vis)
    reveal(vis)
end

end
