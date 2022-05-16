module dlcap
using Parameters
using VoronoiFVM,ExtendableGrids,GridVisualize
using LiquidElectrolytes
using CompositeStructs
using PyPlot,Colors


@siunits nm cm μF

@composite @kwdef mutable struct HalfCellData
    ElectrolyteData...
    TwoElectrodeCell...
end

function main(;n=100,ϕmax=0.9,kwargs...)
    hmin=1.0e-3*nm
    hmax=1.0*nm
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)


    function bcondition(f,u,bnode,data)
        @unpack iϕ,Γ_we,ϕ_we = data
        boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_we,value=ϕ_we)
        bulkbc(f,u,bnode,data)
    end


    celldata=HalfCellData()
    
    cell=NPPSystem(grid;bcondition,celldata)


    molarities=[0.001,0.01,0.1,1]

    vis=GridVisualizer(resolution=(500,300),legend=:rt,clear=true,xlabel="φ/V",ylabel="C_dl/(μF/cm^2)",Plotter=PyPlot)
    hmol=1/length(molarities)
    for imol=1:length(molarities)
	c=RGB(1-imol*hmol,0,imol*hmol)
	t=@elapsed volts,caps=doublelayercap(cell;ϕmax,n,molarity=molarities[imol],kwargs...)
	cdl0=Cdl0(cell.physics.data)
	@info "elapsed=$(t)"
	scalarplot!(vis,volts,caps/(μF/cm^2),
		    color=c,clear=false,label="$(molarities[imol])M",markershape=:none)
	scalarplot!(vis,[0],[cdl0]/(μF/cm^2),clear=false,markershape=:circle,label="")
    end
#    save(plotsdir("1DResults.pdf"),vis)
    reveal(vis)
end

end
