function EquilibriumData(electrolyte::AbstractElectrolyteData)
    EquilibriumData(;
        N=electrolyte.nc,
        T=electrolyte.T,
        p_ref=electrolyte.p_bulk,
        pscale=electrolyte.pscale,
        E_ref=electrolyte.ϕ_bulk,
        n0_ref=ufac"Na"/electrolyte.v0,
        χ=electrolyte.ε-1.0,
        z=electrolyte.z,
        κ=electrolyte.κ,
        molarity=ufac"N_A"*electrolyte.c_bulk[1],
        n_E=ufac"N_A"*electrolyte.c_bulk)
end
