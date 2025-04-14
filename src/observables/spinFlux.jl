function calculate_spinFlux(Ivec::MPS, rho::MPS, sites::Vector{<:Index}, fluxNeededVec::Vector{Int})
    # Initialize vector to store current MPOs
    obsCurrVec = Vector{MPO}(undef, length(fluxNeededVec)-1)
    
    # Current measurement loop
    for i in 1:length(fluxNeededVec)-1
        aobs = OpSum()
        aobs += 4.0, "SxL", fluxNeededVec[i], "SyL", fluxNeededVec[i+1]
        aobs += -4.0, "SyL", fluxNeededVec[i], "SxL", fluxNeededVec[i+1]
        obsCurrVec[i] = MPO(aobs, sites)
    end
    
    # Calculate currents
    currents = Vector{ComplexF64}(undef, length(obsCurrVec))
    for j in 1:length(obsCurrVec)
        numerator = inner(Ivec, Apply(obsCurrVec[j], rho))
        denominator = inner(Ivec, rho)
        currents[j] = (numerator/denominator)
    end
    
    return currents
end