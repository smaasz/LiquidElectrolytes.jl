function halfcellbc(f,u,bnode,data)
    @unpack iϕ,ip,nc,bdata=data
    boundary_dirichlet!(f,u,bnode,species=iϕ,region=bdata.Γ_we,value=bdata.ϕ_we)
    boundary_dirichlet!(f,u,bnode,species=iϕ,region=bdata.Γ_bulk,value=bdata.ϕ_bulk)
    boundary_dirichlet!(f,u,bnode,species=ip,region=bdata.Γ_bulk,value=bdata.p_bulk)
    for ic=1:nc
        boundary_dirichlet!(f,u,bnode,species=ic,region=bdata.Γ_bulk,value=data.c_bulk[ic])
    end
end


@with_kw mutable struct HalfCellBData
    Γ_bulk::Int=2
    Γ_we::Int=1
    ϕ_bulk::Float64=0
    ϕ_we::Float64=0
    p_bulk::Float64=0
end


function ideally_polarizing_halfcell(X::AbstractVector;bdata=HalfCellBData(),electrolyte=ElectrolyteData(bdata), kwargs...)
    grid=simplexgrid(X)
    bdata=HalfCellBData()
    @info bdata
    electrolyte=ElectrolyteData(bdata; kwargs...)
    @info electrolyte
    NPPSystem(grid;bcondition=halfcellbc,electrolyte)
end



const ihplus=1
const iso4=2
const io2=3
const R0=1.0e-6
const β=0.5

function oxhalfcellbc(f,u,bnode,data)
    @unpack iϕ,ip,nc,bdata=data
    halfcellbc(f,u,bnode,data)
    c0,barc=c0_barc(u,data)
    if bnode.region==bdata.Γ_we
        μhplus=xlog(u[ihplus]/barc)#/(R*data.T)
        μo2=xlog(u[io2]/barc)#/(R*data.T)
        r=R0*(exp(-β*(μo2-μhplus)) - exp((1-β)*(μo2-μhplus)))
        f[ihplus]+=r
        f[io2]-=r
    end
end


function oxygen_halfcell(X::AbstractVector;electrolyte=(),kwargs...)
    grid=simplexgrid(X)
    bdata=HalfCellBData()
    @info bdata
    electrolyte=ElectrolyteData(bdata;
                                nc=3,
                                z=[1,-2,0],
                                kwargs...
                                )
    @info electrolyte
    NPPSystem(grid;bcondition=oxhalfcellbc,electrolyte)
end

function voltagesweep(sys;ϕmax=1,δ=1.0e-4,n=100)

    factory=VoronoiFVM.TestFunctionFactory(sys)
    bdata=sys.physics.data.bdata
    
    tf=testfunction(factory,[bdata.Γ_we],[bdata.Γ_bulk] )
    
    dϕ=ϕmax/n
    vplus = zeros(0)
    iplus = zeros(0)
    vminus = zeros(0)
    iminus = zeros(0)

    bdata=boundarydata(sys)
    data=electrolytedata(sys)
    bdata.ϕ_we=0
    

    iϕ=data.iϕ
    inival0 = solve(sys,inival=nppunknowns(sys))
    inival=copy(inival0)
    sol=copy(inival0)
    for dir in [1,-1]
        inival .= inival0
        ϕ = 0.0
        while ϕ <= ϕmax
            @info ϕ
            bdata.ϕ_we=dir * ϕ
            solve!(sol, inival, sys)
            inival .= sol
            I=integrate(sys,tf,sol)
            if dir == 1
                push!(vplus, dir * ϕ)
                push!(iplus, I[ihplus])
            else
                push!(vminus, dir * ϕ)
                push!(iminus, I[ihplus])
            end
            ϕ += dϕ
        end
    end
    vcat(reverse(vminus),vplus),vcat(reverse(iminus),iplus) 
end

    



function doublelayercap(sys;ϕmax=1,δ=1.0e-4,n=100,molarity=0.1,kwargs...)

    dϕ=ϕmax/n
    vplus = zeros(0)
    cdlplus = zeros(0)
    vminus = zeros(0)
    cdlminus = zeros(0)
    bdata=boundarydata(sys)
    data=electrolytedata(sys)
    bdata.ϕ_we=0
    
    data.c_bulk.=molarity*mol/dm^3
    
    iϕ=data.iϕ
    inival0 = solve(sys,inival=nppunknowns(sys))
    inival=copy(inival0)
    sol=copy(inival0)
    for dir in [1,-1]
        inival .= inival0
        ϕ = 0.0
        while ϕ <= ϕmax
            @info ϕ 
            bdata.ϕ_we=dir * ϕ
            solve!(sol, inival, sys; control=VoronoiFVM.SolverControl(;kwargs...))
            Q = integrate(sys, sys.physics.reaction, sol)
            bdata.ϕ_we=dir * ϕ+δ
            inival .= sol
            solve!(sol, inival, sys)
            Qδ = integrate(sys, sys.physics.reaction, sol)
            
            cdl = (Qδ[iϕ] - Q[iϕ]) / δ

            if dir == 1
                push!(vplus, dir * ϕ)
                push!(cdlplus, cdl)
            else
                push!(vminus, dir * ϕ)
                push!(cdlminus, cdl)
            end
            ϕ += dϕ
        end
    end
    vcat(reverse(vminus),vplus),vcat(reverse(cdlminus),cdlplus) 
end

    

