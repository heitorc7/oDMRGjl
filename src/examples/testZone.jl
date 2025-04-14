using ITensors
using ITensorMPS
using Plots
using Statistics
include("../sites/twoSpinHalf.jl")
include("../objects/liouvillian.jl")
include("../objects/Ivec.jl")
include("../observables/spinFlux.jl")
include("../observables/magnonDensity.jl")
include("../analyticalBenchmarking/spinFluxXXX.jl")
include("../analyticalBenchmarking/magnonDensityXXX.jl")

# Physical parameters: system size and coupling factor gamma
L = 5
gamma = 0.1
# Create sites
sites = siteinds("TwoSpinHalf", L)

# Physical Parameters: Jx, Jy, Jz, dissipators location, dissipators temperature, coupling vector, magnetic field
params = (
    Jx = 1.0, Jy = 1.0, Jz = 1.0,
    dissipatorsVec = [1, L],
    dissipatorsTempValues = [1.0, 0.0],
    gammaVec = [gamma, gamma],
    hVec = zeros(L)
)

# Construct LdL MPO
LdL = LdLXYZConstruct(sites, params...; maxdim=500, cutoff=1e-8)
rho = make_ivec(L, sites)

# Plan to do 5 passes or 'sweeps' of DMRG,
# setting maximum MPS internal dimensions
# for each sweep and maximum truncation cutoff
# used when adapting internal dimensions:
nsweeps = 7000
# maxdim = [1]
cutoff = 1E-10
maxdim = [10,20,100,100,200, 200, 200, 200]
# maxdim=1

# Run the DMRG algorithm, returning energy
# (dominant eigenvalue) and optimized MPS
# for i in (100)
energy, rho = dmrg(LdL,rho; nsweeps, maxdim, cutoff, outputlevel=0)
println("Final energy = $energy")

# Now we calculate the flux
Ivec = make_ivec(L, sites)
currents = calculate_spinFlux(Ivec, rho, sites, [1:L;])

println("Analytical value of the current: ", -fluxXXX(L, gamma))
println("Numerical approximation of the current: ", mean(currents))

analyticalMag = magnonDensityXXX(L, gamma)
magDensity = calculate_magnetization(Ivec, rho, sites)

println("Analytical Magnetization profile: ", analyticalMag)
println("Numerical magProf Approx (real) : ", real.(magDensity))
println("Numerical magProf Approx (cmpx) : ", imag.(magDensity))
plot(real.(analyticalMag))
plot!(real.(magDensity), linestyle=:dash)
gui()

# Keeps plots open
read(stdin, Char)