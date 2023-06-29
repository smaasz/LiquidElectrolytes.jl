#=
# Fe redox half cell
([source code](SOURCE_URL))

I-V sweep for ``Fe^{2+} \to Fe^{3+} + e^-``

Methods called:
- [`ElectrolyteData`](@@ref)
- [`ivsweep`](@@ref)
- [`dlcapsweep`](@@ref)
- [`PNPSystem`](@@ref)

=#

module Example110_Fe23Cell
using LessUnitful
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors
using StaticArrays
using InteractiveUtils

function main(;nref=0,
              compare=false,
              eneutral::Bool=false,
              voltages=(-1:0.01:1)*ufac"V",
              dlcap=false,
              R0=1.0e-10,
              molarities=[0.001,0.01,0.1,1],
              scheme=:μex,
              xmax=1,
              κ=10.0,
              Plotter=PyPlot,
              new=false,
              kwargs...)

    @local_phconstants N_A e R ε_0
    F=N_A*e
    @local_unitfactors cm μF mol dm s mA A nm


    
    defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose="e",
              reltol=1.0e-8,
              tol_mono=1.0e-10)

    kwargs=merge(defaults, kwargs) 

    hmin=1.0e-1*nm*2.0^(-nref)
    hmax=1.0*nm*2.0^(-nref)
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)


    R0=R0*ufac"mol/(cm^2*s)"
    Δg = 0.0
    β = 0.5
    ihplus=1
    ife2 = 2
    ife3 = 3
    iso4 = 4

    function halfcellbc(f,u,bnode,data)
        (;nc,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,RT)=data
        bulkbcondition(f,u,bnode,data;region=Γ_bulk)
        if bnode.region==Γ_we
            if !data.eneutral
                boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
            end
            c0,barc=c0_barc(u,data)
            μfe2=chemical_potential(u[ife2], barc, u[ip], v[ife2]+κ*v0, data)
	    μfe3=chemical_potential(u[ife3], barc, u[ip], v[ife2]+κ*v0, data)
            r=rrate(R0,β,(μfe2 - μfe3 + Δg - data.eneutral*F*(u[iϕ]-ϕ_we))/RT)
            f[ife2]-=r
            f[ife3]+=r
        end
        nothing
    end
    
    
    celldata=ElectrolyteData(;nc=4,
                             z=[1,2,3,-2],
                             eneutral,
                             κ=fill(κ,4),
                             Γ_we=1,
                             Γ_bulk=2,
                             scheme)

    (;iϕ::Int,ip::Int)=celldata
    
    celldata.c_bulk[ihplus]=1.0*mol/dm^3
    celldata.c_bulk[ife2]=0.1*mol/dm^3
    celldata.c_bulk[ife3]=0.1*mol/dm^3
    celldata.c_bulk[iso4]=0.75*mol/dm^3

    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-12)
    
    cell=PNPSystem(grid;bcondition=halfcellbc,celldata)
    
    ## Compare electroneutral and double layer cases
    if compare

        celldata.eneutral=false
	volts,currs, sols=ivsweep(cell;voltages,ispec=ife2,kwargs...)

        celldata.eneutral=true
        nvolts,ncurrs, nsols=ivsweep(cell;voltages,ispec=ife2,kwargs...)

        vis=GridVisualizer(;Plotter,resolution=(600,400),clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(A/m^2)")
        scalarplot!(vis,volts,-currs,color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP")
        scalarplot!(vis,nvolts,-ncurrs,clear=false,color=:green,markershape=:none,label="NNP")
        return reveal(vis)
    end
    

    ## Calculate double layer capacitances
    if dlcap
        
        vis=GridVisualizer(;Plotter,resolution=(500,300),legend=:rt,clear=true,xlabel="φ/V",ylabel="C_dl/(μF/cm^2)")
        hmol=1/length(molarities)
        for imol=1:length(molarities)
            color=RGB(1-imol/length(molarities),0,imol/length(molarities))
	    volts,caps=dlcapsweep(cell;voltages,molarity=molarities[imol],kwargs...)
	    scalarplot!(vis,volts,caps/(μF/cm^2);
                        color,
		        clear=false,
                        label="$(molarities[imol])M")
        end
        return  reveal(vis)
    end
    
    ## Full calculation

    volts,currs, sols=ivsweep(cell;voltages,ispec=ife2,kwargs...)
    tsol=VoronoiFVM.TransientSolution(sols,volts)

    for it=1:length(tsol.t)
        tsol.u[it][ife2,:]/=mol/dm^3
        tsol.u[it][ife3,:]/=mol/dm^3
    end

    xmax=xmax*nm
    xlimits=[0,xmax]
    vis=GridVisualizer(;Plotter,resolution=(1200,400),layout=(1,5),clear=true)
    aspect=3.5*xmax/(tsol.t[end]-tsol.t[begin])
    scalarplot!(vis[1,1],F*currs/(mA/cm^2),volts,markershape=:none,title="IV",xlabel="I",ylabel="ϕ")
    scalarplot!(vis[1,2],cell,tsol;species=ife2,aspect,xlimits,title="Fe2+",colormap=:summer,ylabel="ϕ")
    scalarplot!(vis[1,3],cell,tsol;species=ife3,aspect,xlimits,title="Fe3+",colormap=:summer,ylabel="ϕ")
    scalarplot!(vis[1,4],cell,tsol;species=iϕ,aspect,xlimits,title="ϕ",colormap=:bwr,ylabel="ϕ")
    scalarplot!(vis[1,5],cell,tsol;species=ip,aspect,xlimits,title="p",colormap=:summer,ylabel="ϕ")

    reveal(vis)
end

end

#=
```@example Example110_Fe23Cell
Example110_Fe23Cell.main()
```
=#


#=
```@example Example110_Fe23Cell
Example110_Fe23Cell.main(compare=true)
```
=#
