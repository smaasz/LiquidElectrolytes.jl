function EquilibriumData(electrolyte::AbstractElectrolyteData)
    EquilibriumData(
    N=electrolyte.nc,
    T=electrolyte.T,
    p_ref=electrolyte.p_bulk,
    pscale=electrolyte.pscale,
    E_ref=electrolyte.ϕ_bulk,
    n0_ref=N_A/electrolyte.v0,
    χ=electrolyte.ε-1.0,
    z=electrolyte.z,
    κ=electrolyte.κ,
    molarity=N_A*electrolyte.c_bulk[1],
    n_E=N_A*electrolyte.c_bulk)
end
