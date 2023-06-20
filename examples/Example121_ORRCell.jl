#=
```@meta
Draft=true
```


# ORR cell
([source code](SOURCE_URL))

I-V sweep for Oxygen Redox

O2+4H^+ + 4e^- ---> 2H_2O

Methods called:
- [`ElectrolyteData`](@@ref)
- [`ivsweep`](@@ref)
- [`dlcapsweep`](@@ref)
- [`PNPSystem`](@@ref)

=#

module Example121_ORRCell
using ExtendableGrids, GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using PyPlot, Colors
using StaticArrays
using LessUnitful


function main(;
    voltages = -1:0.05:1,
    compare = false,
    molarity = 0.1,
    nref = 0,
    κ = 10.0,
    eneutral = false,
    scheme = :μex,
    Plotter = PyPlot,
    R0::Float64 = 4.0e-15,
    epsreg = 1.0e-20,
    kwargs...,
)

    @local_phconstants R N_A e
    @local_unitfactors nm cm μF mol dm s
    F = N_A * e


    defaults = (;
        max_round = 3,
        tol_round = 1.0e-10,
        reltol = 1.0e-7,
        tol_mono = 1.0e-7,
        verbose = "e",
    )
    kwargs = merge(defaults, kwargs)

    hmin = 1.0e-1 * nm * 2.0^(-nref)
    hmax = 1.0 * nm * 2.0^(-nref)
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)


    R0 = R0 * mol / (cm^2 * s)
    Δg = 0.0
    β = 0.5
    ϕ_we = 0.0
    ihplus = 1
    iso4 = 2
    io2 = 3
    ih2o = 4
    z = [1, -2, 0, 0]
    κ_r=κ
    κ = [κ, κ, 0, 0]

    function halfcellbc(f, u, bnode, data)
        bulkbcondition(f, u, bnode, data)
        (; iϕ, eneutral, ϕ_we, Γ_we, RT) = data

        if bnode.region == Γ_we
            f .= 0.0
            if !eneutral
                boundary_dirichlet!(
                    f,
                    u,
                    bnode;
                    species = iϕ,
                    region = data.Γ_we,
                    value = data.ϕ_we,
                )
            end
            μsolv, μ = @inline  chemical_potentials!(MVector{5,eltype(u)}(undef), u, data)
            A = (4 * μ[ihplus] + μ[io2] - 2μ[ih2o] + Δg + eneutral * 4 * F * (u[iϕ] - data.ϕ_we)) / (RT)
            r = rrate(R0, β, A)
            f[ihplus] -= 4 * r
            f[io2] -= r
            f[ih2o] += 2r + 4κr # shedding of the solvation shell of ihplus
        end
    end




    celldata =
        ElectrolyteData(; nc = 4, z, κ, Γ_we = 1, Γ_bulk = 2, eneutral, scheme, epsreg)

    (; iϕ, c_bulk) = celldata

    c_slv=1.0e-5* mol / dm^3
    c_bulk[io2] = 0.001 * mol / dm^3
    c_bulk[iso4] = molarity * mol / dm^3
    c_bulk[ihplus] = 2.0 * molarity * mol / dm^3
    c_bulk[ih2o] = 55.4*mol/dm^3 - c_bulk[io2]  - (1+κ[iso4])*c_bulk[iso4] - (1+κ[ihplus])*c_bulk[ihplus] - c_slv


    @assert isapprox(celldata.c_bulk' * celldata.z, 0, atol = 1.0e-12)

    cell = PNPSystem(grid; bcondition = halfcellbc, celldata)
    factory = TestFunctionFactory(cell)
    tf = testfunction(factory, [celldata.Γ_we], [celldata.Γ_bulk])

    ###########################################################################
    ## Compare electroneutral and double layer cases
    if compare
        celldata.eneutral = false
        volts, currs, sols = ivsweep(cell; voltages, ispec = io2, kwargs...)

        celldata.eneutral = true
        nvolts, ncurrs, nsols = ivsweep(cell; voltages, ispec = io2, kwargs...)

        vis = GridVisualizer(;
            Plotter,
            resolution = (600, 400),
            clear = true,
            legend = :lt,
            xlabel = "Δϕ/V",
            ylabel = "I/(A/m^2)",
        )
        scalarplot!(
            vis,
            volts,
            currs,
            color = "red",
            markershape = :utriangle,
            markersize = 7,
            markevery = 10,
            label = "PNP",
        )
        scalarplot!(
            vis,
            nvolts,
            ncurrs,
            clear = false,
            color = :green,
            markershape = :none,
            label = "NNP",
        )
        return reveal(vis)
    end

    ###########################################################################
    vis = GridVisualizer(resolution = (1200, 400), layout = (2, 3), Plotter = Plotter)
    
    volts, currs, sols = ivsweep(cell; voltages, ispec = io2, kwargs...)
    tsol = VoronoiFVM.TransientSolution(sols, volts)
    Io2=zeros(length(tsol.t))
    Ih2o=zeros(length(tsol.t))
    for it = 1:length(tsol.t)
        celldata.ϕ_we=volts[it]
        I = integrate(cell, tf, tsol.u[it])
        Io2[it]=I[io2]
        Ih2o[it]=I[ih2o]
        tsol.u[it][io2, :] /= mol / dm^3
        tsol.u[it][ihplus, :] /= mol / dm^3
        tsol.u[it][iso4, :] /= mol / dm^3
        tsol.u[it][ih2o, :] /= mol / dm^3
    end
#    return Io2, Ih2o, volts, currs
    xmax = 3 * nm
    xlimits = [0, xmax]
    aspect = 3.5 * xmax / (tsol.t[end] - tsol.t[begin])

    scalarplot!(
        vis[1, 1],
        volts,
        -22*Io2,
        markershape = :none,
        label = "Io2",
        xlabel = "V",
        color=:red,
    )
    scalarplot!(
        vis[1, 1],
        volts,
        Ih2o,
        markershape = :none,
        label = "Ih2o",
        xlabel = "V",
        legend = :lt,
        color=:blue,
        clear=false,
    )

    scalarplot!(
        vis[2, 1],
        cell,
        tsol;
        species = io2,
        aspect,
        xlimits,
        title = "O2",
        colormap = :summer,
    )
    scalarplot!(
        vis[2, 2],
        cell,
        tsol;
        species = ihplus,
        aspect,
        xlimits,
        title = "H+",
        colormap = :summer,
    )
    scalarplot!(
        vis[2, 3],
        cell,
        tsol;
        species = ih2o,
        aspect,
        xlimits,
        title = "H20",
        colormap = :summer,
    )
    #    scalarplot!(vis[1,4],cell,tsol;species=iϕ,aspect,xlimits,title="ϕ",colormap=:bwr)
    scalarplot!(
        vis[1, 2],
        tsol[io2, 1, :],
        volts,
        label = "O2",
        xlabel = "c",
        color = :green,
        clear = false,
        legend = :rc,
    )
    scalarplot!(
        vis[1, 3],
        tsol[ihplus, 1, :],
        volts,
        title = "c(0)",
        xlabel = "c",
        ylabel = "V",
        label = "H+",
        color = :red,
    )
    scalarplot!(
        vis[1, 3],
        tsol[iso4, 1, :],
        volts,
        label = "SO4--",
        color = :blue,
        clear = false,
        legend = :rc,
    )

    reveal(vis)
end

end
