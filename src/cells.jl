function doublelayercap(sys;vrange=-1:0.05:1,δ=1.0e-4,molarity=0.1,solver_kwargs...)
    ranges=splitz(vrange)

    allprogress=sum(length,ranges)
    vplus = zeros(0)
    cdlplus = zeros(0)
    vminus = zeros(0)
    cdlminus = zeros(0)
    data=electrolytedata(sys)
    data.ϕ_we=0
    
    data.c_bulk.=molarity*mol/dm^3
    
    iϕ=data.iϕ
    inival0 = solve(sys,inival=pnpunknowns(sys))
    inival=copy(inival0)
    sol=copy(inival0)
    ϕprogress = 0

    control=VoronoiFVM.SolverControl(max_round=3, tol_round=1.0e-10;solver_kwargs...)
    @withprogress for range in ranges
        sol .= inival0
        for ϕ in range
            try
                data.ϕ_we= ϕ
                inival.=sol
                solve!(sol, inival, sys; control)
            catch e
                @warn "δ=0 ϕ_we=$(data.ϕ_we)"
                rethrow(e)
            end

            Q = integrate(sys, sys.physics.reaction, sol)
            
            try
                data.ϕ_we=ϕ+δ
                inival .= sol
                solve!(sol, inival, sys; control)
            catch e
                @warn "δ=δ ϕ_we=$(data.ϕ_we)"
                rethrow(e)
            end
            Qδ = integrate(sys, sys.physics.reaction, sol)
            cdl = (Qδ[iϕ] - Q[iϕ]) / δ

            if range[end]>range[1]
                push!(vplus, ϕ)
                push!(cdlplus, cdl)
            else
                push!(vminus, ϕ)
                push!(cdlminus, cdl)
            end
            ϕprogress +=1
            @logprogress ϕprogress/allprogress
        end
    end
    vcat(reverse(vminus),vplus),vcat(reverse(cdlminus),cdlplus) 
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
    weights=ones(data.nc+2)
    weights[data.ip]=0
    mynorm=u->wnorm(u,weights,Inf)
    myrnorm=u->wnorm(u,weights,1)
    data=electrolytedata(sys)
    data.ϕ_we=0
    control=SolverControl(;solver_kwargs...)
    iϕ=data.iϕ
    inival0 = solve(sys;inival=pnpunknowns(sys), solver_kwargs...)
    inival=copy(inival0)
    sol=copy(inival0)
    ϕprogress=0
    dirs=[-1,1]
    @withprogress  for dir in dirs
        inival .= inival0
        ϕ = 0.0
        while ϕ <= ϕmax
            data.ϕ_we=dir * ϕ
            solve!(sol, inival, sys;mynorm,myrnorm,control) 
            inival .= sol
#            I=integrate(sys,tf,sol)
            I=-integrate(sys,sys.physics.breaction,sol; boundary=true)[:,data.Γ_we]
            if dir == 1
                push!(splus, copy(sol))
                push!(vplus, dir * ϕ)
                push!(iplus, I[ispec]*F/(mA/cm^2))
            else
                push!(sminus, copy(sol))
                push!(vminus, dir * ϕ)
                push!(iminus, I[ispec]*F/(mA/cm^2))
            end
            ϕ += dϕ
            ϕprogress +=dϕ
            @logprogress ϕprogress/(2ϕmax)
        end
    end
    if dirs[1]==1
        return vplus,iplus,splus
    end
    popfirst!(iminus)
    popfirst!(vminus)
    popfirst!(sminus)
    vcat(reverse(vminus),vplus),vcat(reverse(iminus),iplus), vcat(reverse(sminus),splus) 
end

    
function splitz(range)
    if range[1]>=0
        return [0:step(range):range[end]]
    elseif range[end]<=0
        return [0.0:-step(range):range[1]]
    else
        [0:-step(range):range[1],0:step(range):range[end]]
    end
end

