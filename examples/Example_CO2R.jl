module Example_CO2R
using LessUnitful
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors
using StaticArrays
using InteractiveUtils

ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
using PyCall

function main(;nref=0,
              # compare=false,
              # eneutral::Bool=false,
              voltages=(-1.5:0.1:-0.7)*ufac"V",
              dlcap=false,
              # R0=1.0e-10,
              molarities=[0.001,0.01,0.1,1],
              scheme=:μex,
              xmax=1,
              κ=10.0,
              Plotter=PyPlot,
              new=false,
              kwargs...)

    @local_phconstants N_A e R ε_0
    F=N_A*e
    c=299792458 * ufac"m/s"
    @local_unitfactors cm μF mol dm s mA A nm


    
    defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose="e",
              reltol=1.0e-8,
              tol_mono=1.0e-10)

    kwargs=merge(defaults, kwargs) 

    hmin=1.0e-1*ufac"μm"*2.0^(-nref)
    hmax=1.0*ufac"μm"*2.0^(-nref)
    #L=20.0*nm
    L=80.0 * ufac"μm"
    X=geomspace(0,L,hmin,hmax)
    grid=simplexgrid(X)


    # R0=R0*ufac"mol/(cm^2*s)"
    # Δg = 0.0
    # β = 0.5

    sigmas = []
    energies = []

    ϕs = []
    rs = []
    ϕ_curr = 100.0

    T   = 273.15 + 25 * ufac"K"
    pH  = 6.8

    Hcp_CO  = 9.7e-6 * ufac"mol/(m^3 * Pa)"
    Hcp_CO2 = 3.3e-4 * ufac"mol/(m^3 * Pa)"

    E_ads_CO2 = 0.657600203 * ufac"eV"
    frequencies = [136.85, 183.6, 212.95, 250.7, 306.0, 510.55, 562.25, 1176.05, 1889.85] * ufac"h * c / cm / eV"
    py"""
    from ase.thermochemistry import HarmonicThermo

    def get_thermal_correction_adsorbate(T, frequencies):
        thermo = HarmonicThermo(frequencies)
        return thermo.get_helmholtz_energy(T, verbose=False)
    """
    harmonic_adsorbate_correction   = py"get_thermal_correction_adsorbate"(T, frequencies) * ufac"eV"

    C_gap = 20 * ufac"μF/cm^2"
    ϕ_pzc = 0.16 * ufac"V"
    
    ikplus      = 1
    ihco3       = 2
    ico3        = 3
    ico2        = 4
    ico         = 5
    iohminus    = 6
    ihplus      = 7

    function halfcellbc(f,u,bnode,data)
        (;nc,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,RT, ε)=data

        bulkbcondition(f,u,bnode,data;region=Γ_bulk)


        if bnode.region==Γ_we

            #boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
            
            # Robin b.c. for the Poisson equation

            boundary_robin!(f, u, bnode, iϕ, C_gap / ε, C_gap * (ϕ_we - ϕ_pzc) / ε)

            # Flux conditions for CO2 and CO
            prefactor   = 1.0e8
            γ_CO2       = 1.0
            a_CO2       = γ_CO2 * u[ico2] * ufac"dm^3/mol"
            θ_free      = 0.9999
            
            sigma                           = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)
            electrochemical_correction      = [-0.000286600929, 0.0297720125]' * [(sigma / ufac"μA/cm^2")^2 , (sigma / ufac"μA/cm^2")] * ufac"eV"
            
            ΔG      = (E_ads_CO2 + harmonic_adsorbate_correction + electrochemical_correction) * ph"N_A"
            
            r = prefactor * a_CO2 * θ_free * exp(-ΔG / RT - 2.3*pH)
            
            if ϕ_we == ϕ_curr
                pop!(sigmas)
                pop!(energies)
                pop!(rs)
            else
                ϕ_curr = ϕ_we
                push!(ϕs, ϕ_we)
                #println("$(ϕ_we):$(r)")
            end

            push!(sigmas, sigma)
            push!(energies, ΔG / (ufac"eV" * ph"N_A"))
            #push!(rs, r)
            push!(rs, u[ico2] * ufac"dm^3/mol")

            f[ico2] = r
            f[ico]  = -r
            f[iohminus] = -2*r

            # c0,barc=c0_barc(u,data)
            # μfe2=chemical_potential(u[ife2], barc, u[ip], v[ife2]+κ*v0, data)
	        # μfe3=chemical_potential(u[ife3], barc, u[ip], v[ife2]+κ*v0, data)
            # r=rrate(R0,β,(μfe2 - μfe3 + Δg - data.eneutral*F*(u[iϕ]-ϕ_we))/RT)
            # f[ife2]-=r
            # f[ife3]+=r
        end
        nothing
    end
    
    
    celldata=ElectrolyteData(;nc=7,
                             z=[1,-1,-2,0,0,-1,1],
                             D=[1.957e-9, 1.185e-9, 0.923e-9, 1.91e-9, 2.23e-9, 5.273e-9, 9.310e-9] * ufac"m^2/s", # from Ringe paper
                             T=T,
                             eneutral=false,
                             κ=fill(κ,7),
                             Γ_we=1,
                             Γ_bulk=2,
                             scheme)

    (;iϕ::Int,ip::Int)=celldata
    
    celldata.c_bulk[ikplus]         = 0.1 * mol/dm^3
    celldata.c_bulk[ihco3]          = (0.1 - 9.53936e-8) * mol/dm^3
    celldata.c_bulk[ico3]           = 9.53936e-8 * mol/dm^3
    celldata.c_bulk[ico2]           = 0.033 * mol/dm^3
    celldata.c_bulk[ico]            = 0.0 * mol/dm^3
    celldata.c_bulk[iohminus]       = 10^(pH - 14) * mol/dm^3
    celldata.c_bulk[ihplus]         = 10^(-pH) * mol/dm^3

    println(celldata.c_bulk'*celldata.z)
    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-10)
    
    cell=PNPSystem(grid;bcondition=halfcellbc,celldata)
    
    ## Compare electroneutral and double layer cases
    # if compare

    #     celldata.eneutral=false
	# volts,currs, sols=ivsweep(cell;voltages,ispec=ife2,kwargs...)

    #     celldata.eneutral=true
    #     nvolts,ncurrs, nsols=ivsweep(cell;voltages,ispec=ife2,kwargs...)

    #     vis=GridVisualizer(;Plotter,resolution=(600,400),clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(A/m^2)")
    #     scalarplot!(vis,volts,-currs,color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP")
    #     scalarplot!(vis,nvolts,-ncurrs,clear=false,color=:green,markershape=:none,label="NNP")
    #     return reveal(vis)
    # end

    if dlcap ## Calculate double layer capacitances
        
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
    else     ## Calculate current density-voltage curve
        volts, currs, sols = ivsweep(cell; voltages, ispec=iohminus, kwargs...)
        vis=GridVisualizer(;Plotter, layout=(1,3))
        scalarplot!(vis[1,1],volts,currs*ufac"cm^2/mA",color="red",markershape=:utriangle,markersize=7, markevery=10,label="PNP",clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="I/(mA/cm^2)", yscale=:log)
        scalarplot!(vis[1,2], sigmas, energies, color="black",clear=true,xlabel="σ/(μC/cm^s)",ylabel="ΔE/eV")
        scalarplot!(vis[1,3], ϕs, rs, xlimits=(-1.5,-.60), yscale=:log, xlabel="Δϕ/V", ylabel="c(CO2)/M")
        for r in rs
            println(r)
        end
        return reveal(vis)
    end

    # ## Full calculation

    # volts,currs, sols=ivsweep(cell;voltages,ispec=ife2,kwargs...)
    # tsol=VoronoiFVM.TransientSolution(sols,volts)

    # for it=1:length(tsol.t)
    #     tsol.u[it][ife2,:]/=mol/dm^3
    #     tsol.u[it][ife3,:]/=mol/dm^3
    # end

    # xmax=xmax*nm
    # xlimits=[0,xmax]
    # vis=GridVisualizer(;Plotter,resolution=(1200,400),layout=(1,5),clear=true)
    # aspect=3.5*xmax/(tsol.t[end]-tsol.t[begin])
    # scalarplot!(vis[1,1],F*currs/(mA/cm^2),volts,markershape=:none,title="IV",xlabel="I",ylabel="ϕ")
    # scalarplot!(vis[1,2],cell,tsol;species=ife2,aspect,xlimits,title="Fe2+",colormap=:summer,ylabel="ϕ")
    # scalarplot!(vis[1,3],cell,tsol;species=ife3,aspect,xlimits,title="Fe3+",colormap=:summer,ylabel="ϕ")
    # scalarplot!(vis[1,4],cell,tsol;species=iϕ,aspect,xlimits,title="ϕ",colormap=:bwr,ylabel="ϕ")
    # scalarplot!(vis[1,5],cell,tsol;species=ip,aspect,xlimits,title="p",colormap=:summer,ylabel="ϕ")

    # reveal(vis)
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
