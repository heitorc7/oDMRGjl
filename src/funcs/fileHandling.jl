using ITensors
using HDF5
function fetchname(L, Delta, DeltaG, gammao, sweep, bd, iteration = 1, Jx = 1, Jy = 1, Jz = 0.5)
    output_dir = "simulationData/MPSdata/N_$(L)/"
    mkpath(output_dir)  # Creates directory if not exists
    output_file = joinpath(output_dir, "TransistorData_N-$(L)_gamma-$(gammao)_deltaG-$(DeltaG)_Delta-$(Delta)_sweep-$(sweep)_bd-$(bd)_iter-$(iter).h5")
    return output_file
end

function saveMPS(rho, name)

    L = length(rho)
    # output_file = fetchname(L, Delta, DeltaG, gammao, sweep, bd, iteration = 1, Jx = 1, Jy = 1, Jz = 0.5)
    output_dir = "simulationData/MPSdata/N_$(L)/"
    mkpath(output_dir)  # Creates directory if not exists
    output_file = joinpath(output_dir, name*".h5")

    f = h5open(output_file,"w")
    write(f,"psi",rho)
    close(f)

    # params = (
    # Delta,
    # DeltaG,
    # gammao
    # )
    # return 
end