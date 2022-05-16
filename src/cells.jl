function bulkbc(f,u,bnode,data)
    @unpack iϕ,ip,nc,Γ_bulk,ϕ_bulk,p_bulk,c_bulk=data
    if bnode.region==Γ_bulk
        boundary_dirichlet!(f,u,bnode,species=iϕ,region=Γ_bulk,value=ϕ_bulk)
        boundary_dirichlet!(f,u,bnode,species=ip,region=Γ_bulk,value=p_bulk)
        for ic=1:nc
            boundary_dirichlet!(f,u,bnode,species=ic,region=Γ_bulk,value=data.c_bulk[ic])
        end
    end
end




@composite @kwdef mutable struct TwoElectrodeCell
    Γ_we::Int=1
    ϕ_we::Float64=0.0
    Γ_bulk::Int=2
end


function voltagesweep(sys;ϕmax=0.5,ispec=1,n=100,solver_kwargs...)

    factory=VoronoiFVM.TestFunctionFactory(sys)
    data=sys.physics.data
    
    tf=testfunction(factory,[data.Γ_bulk],[data.Γ_we] )

    dϕ=ϕmax/n
    vplus = zeros(0)
    iplus = zeros(0)
    splus = []
    vminus = zeros(0)
    iminus = zeros(0)
    sminus = []

    data=electrolytedata(sys)
    data.ϕ_we=0
    control=SolverControl(;solver_kwargs...)

    iϕ=data.iϕ
    inival0 = solve(sys;inival=nppunknowns(sys), solver_kwargs...)
    inival=copy(inival0)
    sol=copy(inival0)
    for dir in [1,-1]
        inival .= inival0
        ϕ = 0.0
        while ϕ <= ϕmax
            @info "ϕ=$(ϕ)"
            data.ϕ_we=dir * ϕ
            solve!(sol, inival, sys; control) 
            inival .= sol
            I=integrate(sys,tf,sol)
            if dir == 1
                
                push!(splus, copy(sol))
                push!(vplus, dir * ϕ)
                push!(iplus, I[ispec])
            else
                push!(sminus, copy(sol))
                push!(vminus, dir * ϕ)
                push!(iminus, I[ispec])
            end
            ϕ += dϕ
        end
    end
    popfirst!(iminus)
    popfirst!(vminus)
    popfirst!(sminus)
    
    vcat(reverse(vminus),vplus),vcat(reverse(iminus),iplus), vcat(reverse(sminus),splus) 
end

    


function doublelayercap(sys;ϕmax=1,δ=1.0e-4,n=100,molarity=0.1,solver_kwargs...)
    dϕ=ϕmax/n
    vplus = zeros(0)
    cdlplus = zeros(0)
    vminus = zeros(0)
    cdlminus = zeros(0)
    data=electrolytedata(sys)
    data.ϕ_we=0
    
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
            data.ϕ_we=dir * ϕ
            solve!(sol, inival, sys; control=VoronoiFVM.SolverControl(;solver_kwargs...))
            Q = integrate(sys, sys.physics.reaction, sol)
            data.ϕ_we=dir * ϕ+δ
            inival .= sol
            solve!(sol, inival, sys; control=VoronoiFVM.SolverControl(;solver_kwargs...))
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

    

