function warmUp(rho::MPS, LdL::MPO, convergenceThreshold = 0.00001, cutoff = 1E-10)
    println("Initializing warm-up routine")
    global energyFin = floatmax(Float64)
    global warmUpTag = true
    # convergenceThreshold = 0.00001 #1E-10
    global warmUpTracker = 0
    while warmUpTag
        global warmUpTracker += 1
        energyIni = energyFin;
        global energyFin, rho = dmrg(LdL,rho; nsweeps=1, maxdim=1, cutoff=cutoff, outputlevel=0)
    
        # println("Difference between energies: ", (energyIni-energyFin))
        if (abs(energyFin) < abs(energyIni)) && abs(abs(energyFin) - abs(energyIni)) < convergenceThreshold
            global warmUpTag = false
            println("EndedWarmUpAfter ", warmUpTracker, " sweeps")
        end
    end
    return rho
end