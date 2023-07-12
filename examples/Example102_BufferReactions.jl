module Example102_BufferReactions
using LessUnitful
using ExtendableGrids,GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot,Colors 
using StaticArrays
using InteractiveUtils
using ForwardDiff
using ProgressLogging

function pHvsweep(sys; voltages = (-0.5:0.1:0.5) * ufac"V", solver_kwargs...)
    ranges = LiquidElectrolytes.splitz(voltages)
    F = ph"N_A*e"
    data = sys.physics.data

    vminus = []
    vplus = []
    sminus = []
    splus = []

    data = electrolytedata(sys)
    data.ϕ_we = 0
    control = SolverControl(;
        verbose = true,
        handle_exceptions = true,
        Δp_min = 1.0e-3,
        Δp = 1.0e-2,
        Δp_grow = 1.2,
        Δu_opt = 1.0e-2,
        unorm = u -> LiquidElectrolytes.wnorm(u, data.weights, Inf),
        rnorm = u -> LiquidElectrolytes.wnorm(u, data.weights, 1),
        solver_kwargs...,
    )
    iϕ = data.iϕ
    @info "Solving for 0V..."
    inival = solve(sys; inival = pnpunknowns(sys), control)

    allprogress = voltages[end] - voltages[begin]
    ϕprogress = 0
    for range in ranges
        @info "pHV sweep from $(range[1])V to $(range[end])V..."
        dir = range[end] > range[1] ? 1 : -1

        psol = nothing
        @withprogress begin

            function pre(sol, ϕ)
                data.ϕ_we = dir * ϕ
            end

            function post(sol, oldsol, ϕ, Δϕ)
                ϕprogress += abs(Δϕ)
                @logprogress ϕprogress / allprogress
            end

            function delta(sys, u, v, t, Δt)
                n = LiquidElectrolytes.wnorm(u, data.weights, Inf) * data.v0
            end

            psol = solve(
                sys;
                inival,
                embed = dir * range,
                control,
                pre,
                post,
                delta,
                store_all = true,
            )
        end

        if dir == 1
            vplus = psol.t
            splus = psol.u
        else
            vminus = -psol.t
            sminus = psol.u

            popfirst!(vminus)
            popfirst!(sminus)
        end
    end

    volts = vcat(reverse(vminus), vplus)
    sols = vcat(reverse(sminus), splus)
    VoronoiFVM.TransientSolution(sols, volts)
end

function main(;nref=0,
              voltages=(-1.5:0.1:-0.6)*ufac"V",
              scheme=:μex,
              κ=10.0,
              Plotter=PyPlot,
              kwargs...)

    @local_phconstants N_A e R ε_0 k_B c_0
    F = N_A*e
    @local_unitfactors cm μF mol dm s mA A nm bar eV μA μm


    
    defaults=(; max_round=3,
              tol_round=1.0e-9,
              verbose="e",
              reltol=1.0e-8,
              tol_mono=1.0e-10)

    kwargs=merge(defaults, kwargs) 

    hmin    = 1.0e-1*μm*2.0^(-nref)
    hmax    = 1.0*μm*2.0^(-nref)
    L       = 80.0 * μm
    X       = geomspace(0,L,hmin,hmax)
    grid    = simplexgrid(X)

    # environment parameters
    T   = 273.15 + 25 * ufac"K"
    pH  = 6.8

    # Henry constants
    Hcp_CO  = 9.7e-6 * ufac"mol/(m^3 * Pa)"
    Hcp_CO2 = 3.3e-4 * ufac"mol/(m^3 * Pa)"

    # buffer equations
    ## in base
    ## CO2 + OH- <=> HCO3-
    kbe1 = 4.44e7 / (mol/dm^3)
    kbf1 = 5.93e3 / (mol/dm^3) / ufac"s"
    kbr1 = kbf1 / kbe1
    ## HCO3- + OH- <=> CO3-- + H2O
    kbe2 = 4.66e3 / (mol/dm^3)
    kbf2 = 1.0e8 / (mol/dm^3) / ufac"s"
    kbr2 = kbf2 / kbe2

    ## in acid
    ## CO2 + H20 <=> HCO3- + H+
    kae1 = 4.44e-7 * (mol/dm^3)
    kaf1 = 3.7e-2 / ufac"s"
    kar1 = kaf1 / kae1
    ## HCO3- <=> CO3-- + H+ 
    kae2 = 4.66e-5 / (mol/dm^3)
    kaf2 = 59.44e3 / (mol/dm^3) / ufac"s"
    kar2 = kaf2 / kae2
    ## autoprotolyse
    kwe  = 1.0e-14 * (mol/dm^3)^2
    kwf  = 2.4e-5 * (mol/dm^3) / ufac"s"
    kwr  = kwf / kwe


    C_gap = 20 * ufac"μF/cm^2"
    ϕ_pzc = 0.16 * ufac"V"
    
    ikplus      = 1
    ihco3       = 2
    ico3        = 3
    ico2        = 4
    ico         = 5
    iohminus    = 6
    ihplus      = 7

    function halfcellbc(f,u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, bnode,data) where {Tval, Tv, Tc, Tp, Ti}
        (;nc,Γ_we,Γ_bulk,ϕ_we,ip,iϕ,v,v0,T, RT, ε)=data

        bulkbcondition(f,u,bnode,data;region=Γ_bulk)

        if bnode.region==Γ_we

            #boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)
            
            # Robin b.c. for the Poisson equation
            boundary_robin!(f, u, bnode, iϕ, C_gap / ε, C_gap * (ϕ_we - ϕ_pzc) / ε)

            # homogeneous Neumann conditions, no flux
            # for ic in 1:nc
            #     f[ic] = 0
            # end
        end
        nothing
    end

    function reaction(f, u::VoronoiFVM.NodeUnknowns{Tv, Tc, Tp, Ti}, node, data) where {Tv, Tc, Tp, Ti}  
        # buffer reactions
        rates       = zeros(Tv, 5)
        ## in base
        ## CO2 + OH- <=> HCO3-
        rates[1]    = kbf1 * u[ico2] * u[iohminus]  - kbr1 * u[ihco3]  
        ## HCO3- + OH- <=> CO3-- + H2O
        rates[2]    = kbf2 * u[ihco3] * u[iohminus] - kbr2 * u[ico3]

        ## in acid
        ## CO2 + H20 <=> HCO3- + H+
        rates[3]    = kaf1 * u[ico2]  - kar1 * u[ihco3] * u[ihplus]  
        ## HCO3- <=> CO3-- + H+ 
        rates[4]    = kaf2 * u[ihco3] - kar2 * u[ico3] * u[ihplus]  

        ## autoprotolyse
        rates[5]    = kwf - kwr * u[ihplus] * u[iohminus]  

        #println("$(ForwardDiff.value.(rates))")

        f[ihco3]    += rates[1] - rates[2] + rates[3] - rates[4]
        f[ico3]     += rates[2] + rates[4]
        f[ihplus]   += rates[3] + rates[4] + rates[5]
        f[iohminus] += -rates[1] -rates[2] + rates[5]

        #println("$(ForwardDiff.value.(u))")
        #println("$(node.index): $(ForwardDiff.value.(f))")
        nothing
    end
    
    celldata=ElectrolyteData(;nc=7,
                             na=0,
                             z=[1,-1,-2,0,0,-1,1],
                             D=[1.957e-9, 1.185e-9, 0.923e-9, 1.91e-9, 2.23e-9, 5.273e-9, 9.310e-9] * ufac"m^2/s", # from Ringe paper
                             T=T,
                             eneutral=false,
                             κ=fill(κ,7),
                             Γ_we=1,
                             Γ_bulk=2,
                             scheme)

    (;iϕ::Int,ip::Int)=celldata
    
    celldata.c_bulk[ikplus]         = 0.0 * mol/dm^3
    celldata.c_bulk[ihco3]          = 0.091 * mol/dm^3
    celldata.c_bulk[ico3]           = 2.68e-5 * mol/dm^3
    celldata.c_bulk[ico2]           = 0.033 * mol/dm^3
    celldata.c_bulk[ico]            = 0.0 * mol/dm^3
    celldata.c_bulk[iohminus]       = 10^(pH - 14) * mol/dm^3
    celldata.c_bulk[ihplus]         = 10^(-pH) * mol/dm^3

    celldata.c_bulk[ikplus]         = -celldata.c_bulk'*celldata.z

    #println("$(celldata.c_bulk'*celldata.z)")
    @assert isapprox(celldata.c_bulk'*celldata.z,0, atol=1.0e-10)
    
    cell    = PNPSystem(grid; bcondition=halfcellbc, reaction=reaction, celldata)
    #println("$(size(unknowns(cell))): $(ForwardDiff.value.(unknowns(cell)'))")

    tsol    = pHvsweep(cell; voltages=voltages, kwargs...) 
    
    vis     = GridVisualizer(;Plotter, layout=(1,1))

    idx = Int(length(tsol.t) / 2)
    scalarplot!(vis[1,1], cell.grid, -log10.(tsol[ihplus, :, idx] / (mol/dm^3)), color="red",markershape=:utriangle,markersize=7, markevery=10,label="$(tsol.t[idx])",clear=true,legend=:lt,xlabel="Δϕ/V",ylabel="-log c(H+)", xscale=:log, yscale=:log)
    return reveal(vis)
end

end

#=
```@example Example111_CO2RCell
Example111_CO2RCell.main()
```
=#