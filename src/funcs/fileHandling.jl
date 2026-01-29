using ITensors
using HDF5
function saveMPS(rho, MPSname)

    L = length(rho)
    output_dir = "simulationData/MPSdata/N_$(L)/"
    mkpath(output_dir)  # Creates directory if not exists
    output_file = string(joinpath(output_dir, MPSname), ".h5")
    
    f = h5open(output_file,"w")
    write(f,"psi",rho)
    close(f)
end

function loadMPS(L, MPSname)

    # L = length(rho)
    output_dir = "simulationData/MPSdata/N_$(L)/"
    # mkpath(output_dir)  # Creates directory if not exists
    output_file = string(joinpath(output_dir, MPSname), ".h5")
    
    f = h5open(output_file,"r")
    # write(f,"psi",rho)
    rho = read(f,"psi",MPS)
    close(f)

    return rho
end