"""
    splitz(range::AbstractRange)

If range contains zero, split it into two parts, one with values <=0 and one with values >=0.
Otherwise, return the range or its reverse, such that first value always is the one with the smallest absolute value.
"""
function splitz(range::AbstractRange)
    if range[1]>=0
        return [range]
    elseif range[end]<=0
        return [reverse(range)]
    else
        [0:-step(range):range[1],0:step(range):range[end]]
    end
end

"""
    splitz(range::Vector)

Version of [`splitz(range::AbstractRange)`](@ref) for vectors.
"""
function splitz(range::Vector)
    if range[1]>=0
        return [vcat([0.0],range)]
    elseif range[end]<=0
        return [vcat([0.0],reverse(range))]
    else
        for i=1:length(range)
            if range[i]≈0.0
                return [ reverse(range[1:i]), range[i:end]]
            elseif i>1 && range[i-1]<0.0 && range[i]>0.0
                return [ vcat([0.0],reverse(range[1:i-1])), vcat([0.0],range[i:end])]
            end
        end
    end

end


"""
    bulkbcondition(f,u,bnode,electrolyte)

Bulk boundary condition for electrolyte: set potential, pressure and concentrations to bulk values.
"""
function bulkbcondition(f,u,bnode,data;region=data.Γ_bulk)
    (; iϕ,ip,nc,ϕ_bulk,p_bulk,c_bulk) = data
    if bnode.region==region
        boundary_dirichlet!(f,u,bnode;species=iϕ,region,value=ϕ_bulk)
        boundary_dirichlet!(f,u,bnode;species=ip,region,value=p_bulk)
        for ic=1:nc
            boundary_dirichlet!(f,u,bnode;species=ic,region,value=data.c_bulk[ic])
        end
    end
end

"""
           dlcapweep(sys;voltages=(-1:0.1:1)*ufac"V",
                                  δ=1.0e-4,
                                  molarity=0.1*ufac"mol/dm^3",
                                  solver_kwargs...)

Calculate double layer capacitances for voltages given in `voltages`.
Returns vector of voltages and vector of double layer capacitances.

Assumptions:
- Only one double layer in the system - close to working electrode
- 1D domain
"""
function dlcapsweep(sys;voltages=(-1:0.1:1)*ufac"V",δ=1.0e-4,molarity=0.1*ufac"mol/dm^3",solver_kwargs...)
    ranges=splitz(voltages)
    vplus = zeros(0)
    cdlplus = zeros(0)
    vminus = zeros(0)
    cdlminus = zeros(0)
    data=electrolytedata(sys)
    data.ϕ_we=0

    data.c_bulk.=molarity

    iϕ=data.iϕ
    inival=solve(sys,inival=pnpunknowns(sys))
    allprogress=sum(length,ranges)
    ϕprogress = 0
    
    function show_error(u,δ)
        @show cond(inv(Diagonal(sys.matrix))*Matrix(sys.matrix))
        @show u[1,1:5]
        @show u[2,1:5]
        @show u[3,1:5]
        @show u[4,1:5]
        @show chemical_potentials!(zeros(data.nc),u[:,1],data)
        c0,barc= c0_barc(u[:,1],data)
        @show c0/barc, u[1,1]/barc,u[2,1]/barc
        @error "bailing out at δ=$(δ) ϕ_we=$(data.ϕ_we)V, molarity=$(molarity)"
    end

    
    control=VoronoiFVM.SolverControl(max_round=3, tol_round=1.0e-9;solver_kwargs...)
    @withprogress for range in ranges
        sol = inival
        success=true
        for ϕ in range 
            try
                data.ϕ_we= ϕ
                sol=solve(sys;inival=sol, control)
            catch e
                println(e)
                show_error(sol,0)
                success=false
            end

            if !success
                break
            end
            Q = integrate(sys, sys.physics.reaction, sol)
            
            try
                data.ϕ_we=ϕ+δ
                sol=solve(sys;inival=sol, control)
            catch e
                println(e)
                show_error(sol,δ)
                success=false
            end
            if !success
                break
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

    


"""
       ivsweep(sys;
                  voltages=(-0.5:0.1:0.5)*ufac"V",
                                 ispec=1,
                                 solver_kwargs...)

Calculate working electrode current corresponding to rate for species `ispec` for each voltage in `voltages`.
Returns vector of voltages   and vector of currents.
"""
function ivsweep(sys;voltages=(-0.5:0.1:0.5)*ufac"V",ispec=1,solver_kwargs...)
    ranges=splitz(voltages)
    F=ph"N_A*e"
    
    factory=VoronoiFVM.TestFunctionFactory(sys)
    data=sys.physics.data
    
    tf=testfunction(factory,[data.Γ_bulk],[data.Γ_we] )

    iplus = zeros(0)
    iminus = zeros(0)
    vminus=[]
    vplus=[]
    sminus=[]
    splus=[]
    weights=ones(data.nc+2)
    weights[data.ip]=0
    mynorm=u->wnorm(u,weights,Inf)
    myrnorm=u->wnorm(u,weights,1)
    data=electrolytedata(sys)
    data.ϕ_we=0
    control=SolverControl(;verbose=true,
                          handle_exceptions=true,
                          Δp_min=1.0e-6,
                          Δp=1.0e-2,
                          Δp_grow=1.2,
                          Δu_opt=5.0e-3,
                          solver_kwargs...)
    iϕ=data.iϕ
    inival = solve(sys;inival=pnpunknowns(sys), solver_kwargs...)

    allprogress=voltages[end]-voltages[begin]
    ϕprogress=0
    @withprogress for range in ranges
        dir=range[end]>range[1] ? 1 : -1

        function pre(sol,ϕ)
            data.ϕ_we=dir*ϕ
        end
        
        function post(sol,oldsol, ϕ, Δϕ)
            I=-integrate(sys,sys.physics.breaction,sol; boundary=true)[:,data.Γ_we]
            if dir>0
                push!(iplus, I[ispec]*F)
            else
                push!(iminus, I[ispec]*F)
            end
            ϕprogress += abs(Δϕ)
            @logprogress ϕprogress/allprogress
        end
        
        function delta(sys,u,v,t, Δt)
            n=mynorm(u-v)*data.v0
        end
        psol=solve(sys;inival,embed=dir*range,control,pre,post,delta,store_all=true,mynorm=mynorm,myrnorm=myrnorm)

        if dir==1
            vplus=psol.t
            splus=psol.u
        else
            vminus=-psol.t
            sminus=psol.u

            popfirst!(iminus)
            popfirst!(vminus)
            popfirst!(sminus)
        end
    end


    vcat(reverse(vminus),vplus),vcat(reverse(iminus),iplus), vcat(reverse(sminus),splus) 
end

    
