using ITensors
using ITensorMPS
using LinearAlgebra: dot

function LNRN_LiouvXYZConstruct(
    Jx::Float64,
    Jy::Float64,
    Jz::Float64,
    sites::Vector{Index{Int64}},
    DissipatorsLoc::Vector{Int},
    f::Vector{Float64},
    gamma::Vector{Float64},
    h::Vector{Float64}
)
    L = length(sites)/2
    ampo = AutoMPO();
    
    # Input validation
    if length(DissipatorsLoc) != length(f)
        error("Misleading dissipators temperature vector while constructing Liouvillian")
    end
    if length(DissipatorsLoc) != length(gamma)
        error("Misleading dissipators to gamma vector while constructing Liouvillian")
    end
    if L != length(h)
        error("Misleading magnetic field vector while constructing Liouvillian")
    end

    println("\nConstructing Liouvillian for a spin chain of size $L using $(length(DissipatorsLoc)) dissipators")
    
    # Construct dissipators
    for (j, site) in enumerate(DissipatorsLoc)
        # D_m interactions
        ampo += gamma[j]*f[j],      "S-", 2*site-1, "S-", 2*site
        ampo += -gamma[j]*f[j]*0.5, "S+", 2*site-1, "S-", 2*site-1
        ampo += -gamma[j]*f[j]*0.5, "S+", 2*site, "S-", 2*site
        
        # D_p interactions
        ampo += gamma[j]*(1.0-f[j]),      "S+", 2*site-1, "S+", 2*site
        ampo += -gamma[j]*(1.0-f[j])*0.5, "S-", 2*site-1, "S+", 2*site-1
        ampo += -gamma[j]*(1.0-f[j])*0.5, "S-", 2*site, "S+", 2*site
    end

    # Construct local Hamiltonian terms
    for j in 1:L
        ampo += -2im*h[j], "Sz", 2*j-1
        ampo += 2im*h[j],  "Sz", 2*j
    end

    # Construct nearest-neighbor interactions
    for j in 1:L-1
        ampo += -4im*Jx, "Sx", 2*j-1, "Sx", 2*j+1
        ampo += 4im*Jx,  "Sx", 2*j,   "Sx", 2*j+2

        ampo += -4im*Jy, "Sy", 2*j-1, "Sy", 2*j+1
        ampo += 4im*Jy,  "Sy", 2*j,   "Sy", 2*j+2

        ampo += -4im*Jz, "Sz", 2*j-1, "Sz", 2*j+1
        ampo += 4im*Jz,  "Sz", 2*j,   "Sz", 2*j+2
    end

    return ampo
end

function LNRN_LiouvXYZDConstruct(
    Jx::Float64,
    Jy::Float64,
    Jz::Float64,
    sites::Vector{Index{Int64}},
    DissipatorsLoc::Vector{Int},
    f::Vector{Float64},
    gamma::Vector{Float64},
    h::Vector{Float64}
)
    L = length(sites)
    ampo = AutoMPO();
    
    # Input validation
    if length(DissipatorsLoc) != length(f)
        error("Misleading dissipators temperature vector while constructing Liouvillian")
    end
    if length(DissipatorsLoc) != length(gamma)
        error("Misleading dissipators to gamma vector while constructing Liouvillian")
    end
    if L != length(h)
        error("Misleading magnetic field vector while constructing Liouvillian")
    end

    println("\nConstructing Liouvillian for a spin chain of size $L using $(length(DissipatorsLoc)) dissipators")
    
    # Construct dissipators
    for (j, site) in enumerate(DissipatorsLoc)
        # D_m interaction
        ampo += gamma[j]*f[j],      "S+", 2*site-1, "S+", 2*site
        ampo += -gamma[j]*f[j]*0.5, "S+", 2*site-1, "S-", 2*site-1
        ampo += -gamma[j]*f[j]*0.5, "S+", 2*site, "S-", 2*site
    
        # D_p interactions
        ampo += gamma[j]*(1.0-f[j]),      "S-", 2*site-1, "S-", 2*site
        ampo += -gamma[j]*(1.0-f[j])*0.5, "S-", 2*site-1, "S+", 2*site-1
        ampo += -gamma[j]*(1.0-f[j])*0.5, "S-", 2*site, "S+", 2*site
    end

    # Construct local Hamiltonian terms
    for j in 1:L
        ampo += 2im*h[j],  "Sz", 2*j-1
        ampo += -2im*h[j], "Sz", 2*j
    end

    # Construct nearest-neighbor interactions
    for j in 1:L-1
        ampo += 4im*Jx,  "Sx", 2*j-1, "Sx", 2*j+1
        ampo += -4im*Jx, "Sx", 2*j,   "Sx", 2*j+2

        ampo += 4im*Jy,  "Sy", 2*j-1, "Sy", 2*j+1
        ampo += -4im*Jy, "Sy", 2*j,   "Sy", 2*j+2

        ampo += 4im*Jz, "Sz",   2*j-1, "Sz", 2*j+1;
        ampo += -4im*Jz, "SzR", 2*j,   "Sz", 2*j+2;
    end

    return ampo
end

function LNRN_LdLXYZConstruct(
    sites::Vector{Index{Int64}},
    Jx::Float64,
    Jy::Float64,
    Jz::Float64,
    dissipatorsVec::Vector{Int},
    dissipatorsTempValues::Vector{Float64},
    gammaVec::Vector{Float64},
    hVec::Vector{Float64};
    maxdim::Int=50000,
    cutoff::Float64=1e-8
)
    # L = length(sites)/2
    
    # Construct the Liouvillian and its dagger
    LiouvAMPO = LNRN_LiouvXYZConstruct(Jx, Jy, Jz, sites, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec)
    LiouvDAMPO = LNRN_LiouvXYZDConstruct(Jx, Jy, Jz, sites, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec)
    
    # Convert to MPOs using site indices
    Liouv = MPO(LiouvAMPO, sites)
    LiouvD = MPO(LiouvDAMPO, sites)
    
    # Construct Ld*L with truncation parameters
    # LiouvD_primed = prime(LiouvD, 3)
    # LdL = contract(Liouv, LiouvD_primed; maxdim, cutoff)
    LdL = prime(LiouvD)*Liouv
    
    # Un-prime the result
    LdL = replaceprime(LdL, 2=>1)
    
    return LdL
end