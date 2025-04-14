function calculate_magnetization(Ivec::MPS, rho::MPS, sites::Vector{<:Index})
    # Initialize vector to store current MPOs
    
    magNeededVec = [1:length(sites);];
    obsMagVec = Vector{MPO}(undef, length(magNeededVec));
    
    # Magnon Density Operator creator loop
    for i in 1:length(magNeededVec)
        aobs = OpSum()
        aobs += 2.0, "SzL", magNeededVec[i];
        # aobs += 2.0, "SzR", magNeededVec[i];
        obsMagVec[i] = MPO(aobs, sites)
    end
    
    # Calculate magnetizations
    magDensity = Vector{ComplexF64}(undef, length(obsMagVec))
    for j in 1:length(obsMagVec)
        numerator = inner(Ivec, Apply(obsMagVec[j], rho))
        denominator = inner(Ivec, rho)
        magDensity[j] = (numerator/denominator)
    end

    return magDensity
end