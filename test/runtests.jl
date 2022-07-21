using Test
using LiquidElectrolytes
using LinearAlgebra
using CompositeStructs,Parameters
using Base: @kwdef
using ExtendableGrids,VoronoiFVM
using LessUnitful

@phconstants N_A
@unitfactors dm nm mol

@testset "cdl0" begin
    ely=ElectrolyteData(c_bulk=fill(0.01*mol/dm^3,2).|>unitfactor)
    @test Cdl0(ely)≈ 0.22846691848825248
    edata=EquilibriumData()
    LiquidElectrolytes.set_molarity!(edata,0.01)
    edata.χ=78.49-1
    @test Cdl0(edata)|>unitfactor≈ 0.22846691848825248
    @test Cdl0(EquilibriumData(ely))≈ 0.22846691848825248
end


@testset "dlcap" begin
    
    @composite @kwdef mutable struct HalfCellData<:AbstractElectrolyteData
        ElectrolyteData...
        Γ_we::Int=1
        ϕ_we::Float64=0.0
    end
    
    function bcondition(f,u,bnode,data::HalfCellData)
        @unpack iϕ,Γ_we,ϕ_we = data
        boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_we,value=ϕ_we)
        bulkbc(f,u,bnode,data)
    end


    hmin=1.0e-1*nm
    hmax=1.0*nm
    L=20.0*nm
    X=geomspace(0,L,hmin,hmax)
    molarity=0.1*ufac"M"
    voltages=collect(-1:0.01:1)*ufac"V"
    δ=1.0e-4
    
    grid=simplexgrid(X)
    κ=[0,0]
    acelldata=HalfCellData(;Γ_we=1, Γ_bulk=2,  scheme=:act,κ)
    acell=PNPSystem(grid;bcondition,celldata=acelldata)
    avolts,acaps=doublelayercap(acell;voltages,molarity,δ)

    μcelldata=HalfCellData(;Γ_we=1, Γ_bulk=2, scheme=:μex,κ)
    μcell=PNPSystem(grid;bcondition,celldata=μcelldata)
    μvolts,μcaps=doublelayercap(μcell;voltages,molarity,δ)
    
    ccelldata=HalfCellData(;Γ_we=1, Γ_bulk=2, scheme=:cent,κ)
    ccell=PNPSystem(grid;bcondition,celldata=ccelldata)
    cvolts,ccaps=doublelayercap(ccell;voltages,molarity,δ)
    
    ecell=create_equilibrium_system(grid,EquilibriumData(acelldata))
    evolts,ecaps=calc_Cdl(ecell,vmax=1.0,molarity=0.1,δV=1.0e-4,nsteps=101)

    @show norm((acaps-ecaps)./ecaps)
    @show norm((acaps-μcaps)./acaps)
    @show norm((acaps-ccaps)./acaps)

    @test Cdl0(acelldata) ≈ Cdl0(EquilibriumData(acelldata))
    @test isapprox(acaps,ecaps,rtol=1.0e-3)
    @test isapprox(acaps,μcaps,rtol=1.0e-10)
    @test isapprox(acaps,ccaps,rtol=1.0e-10)
end    
