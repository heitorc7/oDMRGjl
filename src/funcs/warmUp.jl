function warmUp(rho::MPS, LdL::MPO, convergenceThreshold = 0.00001, cutoff = 1E-10)
    println("Initializing warm-up routine")
    energyFin = floatmax(Float64)
    warmUpTracker = 0
    while true
        warmUpTracker += 1
        energyIni = energyFin;
        energyFin, rho = dmrg(LdL,rho; nsweeps=1, maxdim=1, cutoff=cutoff, outputlevel=0)
    
        # println("Difference between energies: ", (energyIni-energyFin))
        if (abs(energyFin) < abs(energyIni)) && abs(abs(energyFin) - abs(energyIni)) < convergenceThreshold
            println("EndedWarmUpAfter ", warmUpTracker, " sweeps")
            break
        end
    end
    return rho
end