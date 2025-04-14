using ITensors
using ITensorMPS
using LinearAlgebra: dot

function LiouvXYZConstruct(
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
        # D_m interactions
        ampo += gamma[j]*f[j], "SmL", site, "SpR", site
        ampo += -gamma[j]*f[j]*0.5, "SpSmL", site
        ampo += -gamma[j]*f[j]*0.5, "SpSmR", site
        
        # D_p interactions
        ampo += gamma[j]*(1.0-f[j]), "SpL", site, "SmR", site
        ampo += -gamma[j]*(1.0-f[j])*0.5, "SmSpL", site
        ampo += -gamma[j]*(1.0-f[j])*0.5, "SmSpR", site
    end

    # Construct local Hamiltonian terms
    for j in 1:L
        ampo += -2im*h[j], "SzL", j
        ampo += 2im*h[j], "SzR", j
    end

    # Construct nearest-neighbor interactions
    for j in 1:L-1
        ampo += -4im*Jx, "SxL", j, "SxL", j+1
        ampo += 4im*Jx, "SxR", j, "SxR", j+1

        ampo += -4im*Jy, "SyL", j, "SyL", j+1
        ampo += 4im*Jy, "SyR", j, "SyR", j+1

        ampo += -4im*Jz, "SzL", j, "SzL", j+1
        ampo += 4im*Jz, "SzR", j, "SzR", j+1
    end

    return ampo
end

function LiouvXYZDConstruct(
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
        # D_m interactions
        ampo += gamma[j]*f[j], "SmL", site, "SpR", site
        ampo += -gamma[j]*f[j]*0.5, "SpSmL", site
        ampo += -gamma[j]*f[j]*0.5, "SpSmR", site
        
        # D_p interactions
        ampo += gamma[j]*(1.0-f[j]), "SpL", site, "SmR", site
        ampo += -gamma[j]*(1.0-f[j])*0.5, "SmSpL", site
        ampo += -gamma[j]*(1.0-f[j])*0.5, "SmSpR", site
    end

    # Construct local Hamiltonian terms
    for j in 1:L
        ampo += 2im*h[j], "SzL", j
        ampo += -2im*h[j], "SzR", j
    end

    # Construct nearest-neighbor interactions
    for j in 1:L-1
        ampo += 4im*Jx, "SxL", j, "SxL", j+1
        ampo += -4im*Jx, "SxR", j, "SxR", j+1

        ampo += 4im*Jy, "SyL", j, "SyL", j+1
        ampo += -4im*Jy, "SyR", j, "SyR", j+1

        ampo += 4im*Jz, "SzL", j, "SzL", j+1
        ampo += -4im*Jz, "SzR", j, "SzR", j+1
    end

    return ampo
end

function LdLXYZConstruct(
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
    # L = length(sites)
    
    # Construct the Liouvillian and its dagger
    LiouvAMPO = LiouvXYZConstruct(Jx, Jy, Jz, sites, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec)
    LiouvDAMPO = LiouvXYZDConstruct(Jx, Jy, Jz, sites, dissipatorsVec, dissipatorsTempValues, gammaVec, hVec)
    
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