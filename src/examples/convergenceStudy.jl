using ITensors
using ITensorMPS
using Plots
using Statistics
include("../sites/twoSpinHalf.jl")
include("../objects/liouvillian.jl")
include("../objects/Ivec.jl")
include("../observables/spinFlux.jl")
include("../observables/magnonDensity.jl")
include("../funcs/warmUp.jl")
include("../analyticalBenchmarking/spinFluxXXX.jl")
include("../analyticalBenchmarking/magnonDensityXXX.jl")

# Physical parameters: system size and coupling factor gamma
L = 10
gamma = 0.1
# Create sites using customized TwoSpinHalf sites
sites = siteinds("TwoSpinHalf", L)

# Physical Parameters: Jx, Jy, Jz, dissipators location, dissipators temperature, coupling vector, magnetic field
params = (
    Jx = 1.0, Jy = 1.0, Jz = 1.0,
    dissipatorsVec = [1, L],
    dissipatorsTempValues = [1.0, 0.0],
    gammaVec = [gamma, gamma],
    hVec = zeros(L)
)

# LdL MPO constructor
LdL = LdLXYZConstruct(sites, params...; maxdim=500, cutoff=1e-8)
rho = make_ivec(L, sites)

cutoff = 1E-10

warmUp(rho, LdL, 0.0000001)
Ivec = make_ivec(L, sites)

# Now initiating iterative maxDim increase oDMRG
global maxdim2 = 2;
maxDimInc = 2;
sweepBDChangeThresholdValue = 0.001 # good values: 0.0005 ~ 0.001
finalSweep = 50; # usually 1000~2000 works for N = 10. Experiment with more for larger sizes

global energyFin = floatmax(Float64)
for iniSweep in 1:finalSweep
    energyIni = energyFin
    global energyFin, rho = dmrg(LdL,rho; nsweeps=1, maxdim=maxdim2, cutoff=cutoff)

    # Now we calculate the flux
    # currents = calculate_spinFlux(Ivec, rho, sites, [1:L;])
    if (abs(energyFin) < abs(energyIni)) && abs(abs(energyFin) - abs(energyIni))/energyIni < sweepBDChangeThresholdValue
        global maxdim2 += maxDimInc
    end
end
println("Final energy = $energyFin")


currents = calculate_spinFlux(Ivec, rho, sites, [1:L;])
println("Analytical value of the current: ", -fluxXXX(L, gamma))
println("Numerical approximation of the current: ", mean(currents))

analyticalMag = magnonDensityXXX(L, gamma)
magDensity = calculate_magnetization(Ivec, rho, sites)

println("Analytical Magnetization profile: ", analyticalMag)
println("Numerical magProf Approx (real) : ", real.(magDensity))
println("Numerical magProf Approx (cmpx) : ", imag.(magDensity))
plot([1:L;], real.(analyticalMag))
# Using the imaginary part of the magnon density as  uncertainty (clearly an overestimatation)
plot!([1:L;], real.(magDensity),yerr=imag.(magDensity), linestyle=:dash)
gui()

# Keeps plots open
read(stdin, Char)