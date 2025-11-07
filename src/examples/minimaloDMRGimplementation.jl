using ITensors
# import ITensors: permute  # Explicitly import permute from ITensors
using ITensorMPS
# using TimeEvoMPS
using Plots
using Statistics
include("../sites/twoSpinHalf.jl")
include("../objects/liouvillian.jl")
include("../objects/Ivec.jl")
include("../observables/spinFlux.jl")
include("../observables/magnonDensity.jl")
include("../funcs/warmUp.jl")
include("../funcs/trotterTimeEvol.jl")
include("../analyticalBenchmarking/spinFluxXXX.jl")
include("../analyticalBenchmarking/magnonDensityXXX.jl")

# Physical parameters: system size and coupling factor gamma
L = 6
gamma = 1.0
# Create sites using customized TwoSpinHalf sites
sites = siteinds("TwoSpinHalf", L)

# Physical Parameters: Jx, Jy, Jz, dissipators location, dissipators temperature, coupling vector, magnetic field
params = (
    Jx = 1.0, Jy = 1.0, Jz = 1.0,
    dissipatorsVec = [1, L],
    dissipatorsTempValues = [1.0, 0.0],
    gammaVec = [gamma, gamma],
    hVec = Vector{Float64}(1000000*zeros(L))
)

global cutoff = 1E-8
global fluxVec = Float64[]

# LdL MPO constructor
LdL = LdLXYZConstruct(sites, params...; maxdim=10000000, cutoff=cutoff)
rho = make_ivec(L, sites)
Ivec = make_ivec(L, sites)

###### Warming up and DMRG routine

warmUp(rho, LdL, 0.001)

# Now initiating iterative maxDim increase oDMRG
global maxdim2 = 2;
maxDimInc = 2;
sweepBDChangeThresholdValue = 0.001 # good values: 0.0005 ~ 0.001
finalSweep = 1000; # usually 1000~2000 works for N = 10. Experiment with more for larger sizes

global energyFin = floatmax(Float64)
for iniSweep in 1:finalSweep
    energyIni = energyFin
    global energyFin, rho = dmrg(LdL,rho; nsweeps=1, maxdim=maxdim2, cutoff=cutoff)

    # Now we calculate the flux
    # flux = mean(calculate_spinFlux(Ivec, rho, sites, [1:L;]))
    # push!(fluxVec, real.(flux))
    if (abs(energyFin) < abs(energyIni)) && abs(abs(energyFin) - abs(energyIni))/energyIni < sweepBDChangeThresholdValue
        global maxdim2 += maxDimInc
    end
end
println("Final energy = $energyFin")

########################## Directly trying to apply the MPO to the MPS #######################################

tau = 0.1
ttotal = 50

# LiouvAMPO = LiouvXYZConstruct(params.Jx, params.Jy, params.Jz, sites, params.dissipatorsVec, params.dissipatorsTempValues, params.gammaVec, params.hVec)
# Liouv = MPO(LiouvAMPO, sites)
# expiH = exp(tau/2 * LiouvAMPO);
# println(expiH)
for t in 0.0:tau:ttotal
    # Sz = expect(rho, "Sz"; sites=c)
    flux = mean(calculate_spinFlux(Ivec, rho, sites, [1:L;]))
    println("$t $flux")
    push!(fluxVec, real.(flux))
    
    t≈ttotal && break
    
    # global rho = apply(expiH, rho)
    # normalize!(rho)
    # tebd!(rho,expiH,tau,ttotal, maxdim=maxdim)
    # tdvp(Liouv, -1.0im, rho; nsteps=100, kwargs...)
    global rho = tdvp(LdL, -tau, rho; 
    nsteps=1,           # Number of time steps
    # time_step=0.01,      # Size of each time step (dt)
    cutoff=1e-10,        # SVD cutoff
    maxdim=1000,          # Maximum bond dimension
    outputlevel=1        # Verbosity (0=silent, 1=some output, 2=verbose)
)
end


################################## Now adding Trotter time evolution #######################################

# tau = 0.01
# ttotal = 50

# gates = ITensor[]
# for j in 1:(L - 1)
#     s1 = sites[j]
#     s2 = sites[j + 1]

#     ## Ok best way to do this is: create a bathTerm (D_m + D_p), 
#     ## then local Hamiltonian, then nearest neighbor interaction
#     # println(j)
#     # for (bathind, site) in enumerate(params.dissipatorsVec)
#     #     println(bathind)
#     #     println(site)
#     if j in params.dissipatorsVec
#         println("There is a dissipator at site $j")     # For debugging
#         global bathTerm =
#         # D_m interaction
#         # params.gammaVec[j]*params.dissipatorsTempValues[j] * op("SmL", s1) * op("SpR", s2)
#         params.gammaVec[j]*params.dissipatorsTempValues[j] * op("SmLSpR", s1) * op("Id", s2)
#         + (-params.gammaVec[j]*params.dissipatorsTempValues[j]*0.5) * op("SpSmL", s1) * op("Id", s2)
#         + (-params.gammaVec[j]*params.dissipatorsTempValues[j]*0.5) * op("SpSmR", s1) * op("Id", s2)
#         # + (-params.gammaVec[j]*params.dissipatorsTempValues[j]*0.5) * op("Id", s1) * op("SpSmR", s2) # missing term?

#         # D_p interaction
#         # + (params.gammaVec[j]*(1.0-params.dissipatorsTempValues[j])) * op("SpL", s1) * op("SmR", s2)
#         + (params.gammaVec[j]*(1.0-params.dissipatorsTempValues[j])) * op("SpLSmR", s1) * op("Id", s2)
#         + (-params.gammaVec[j]*(1.0-params.dissipatorsTempValues[j])*0.5) * op("SmSpL", s1) * op("Id", s2)
#         + (-params.gammaVec[j]*(1.0-params.dissipatorsTempValues[j])*0.5) * op("SmSpR", s1) * op("Id", s2)
#         # + (-params.gammaVec[j]*(1.0-params.dissipatorsTempValues[j])*0.5) * op("Id", s1) * op("SmSpR", s2) # missing term?
#         # Gj = exp(tau / 2 * bathTerm)
#         # push!(gates, Gj)
#     else
#         println("There is NOT a dissipator at site $j")    # For debugging
#         global bathTerm = op("Id", s1) * op("Id", s2)
#     end
#     # end
#     # print(bathTerm)

#     ## Local Hamiltonian
#     localHam = (-2im*params.hVec[j]) * op("SzL", s1) * op("Id", s2)
#     + (2im*params.hVec[j]) * op("SzR", s1) * op("Id", s2)
#     # Gj = exp(tau / 2 * localHam)
#     # push!(gates, Gj)

#     ## Nearest Neighbor interactions
#     nearestNeighbor = 
#     (-4im*params.Jx) * op("SxL", s1) * op("SxL", s2)
#     +(4im*params.Jx) * op("SxR", s1) * op("SxR", s2)

#     +(-4im*params.Jy) * op("SyL", s1) * op("SyL", s2)
#     +(4im*params.Jy) * op("SyR", s1) * op("SyR", s2)

#     +(-4im*params.Jz) * op("SzL", s1) * op("SzL", s2)
#     +(4im*params.Jz) * op("SzR", s1) * op("SzR", s2)

#     ## Sum everything up
#     hj = bathTerm + localHam + nearestNeighbor
#     # hj = localHam
# 	# Gj = exp(-im * tau / 2 * hj)
#     # Gj = exp(tau / 2 * hj)
#     Gj = exp(tau * hj)
#     # Gj = exp(tau / 2 * nearestNeighbor)
# 	push!(gates, Gj)
# end
# ## Last site
# s1 = sites[L - 1]
# s2 = sites[L]

# if L in params.dissipatorsVec
#     println("There is a dissipator at site $L")
#     global bathTerm =
#     # D_m interaction
#     last(params.gammaVec)*last(params.dissipatorsTempValues) * op("Id", s1) * op("SmLSpR", s2)
#     + (-last(params.gammaVec)*last(params.dissipatorsTempValues)*0.5) * op("Id", s1) * op("SpSmL", s2)
#     + (-last(params.gammaVec)*last(params.dissipatorsTempValues)*0.5) * op("Id", s1) * op("SpSmR", s2) * 

#     # D_p interaction
#     + (last(params.gammaVec)*(1.0-last(params.dissipatorsTempValues))) * op("Id", s1) * op("SpLSmR", s2)
#     + (-last(params.gammaVec)*(1.0-last(params.dissipatorsTempValues))*0.5) * op("Id", s1) * op("SmSpL", s2)
#     + (-last(params.gammaVec)*(1.0-last(params.dissipatorsTempValues))*0.5) * op("Id", s1) * op("SmSpR", s2)
#     # println(bathTerm)
#     # Gj = exp(tau / 2 * bathTerm)
#     # push!(gates, Gj)
# else
#     # println("There is NOT a dissipator at site $L")
#     global bathTerm = op("Id", s1) * op("Id", s2)
# end

# ## Local Hamiltonian
# localHam = (-2im*params.hVec[L]) * op("Id", s1) * op("SzL", s2)
# + (2im*params.hVec[L]) * op("Id", s1) * op("SzR", s2)
# # Gj = exp(tau / 2 * localHam)
# # push!(gates, Gj)

# hj = bathTerm + localHam
# # hj = localHam
# # Gj = exp(-im * tau / 2 * hj)  # Real Time Evolution
# Gj = exp(-tau * hj)         # Imaginary Time Evolution
# push!(gates, Gj)

# # append!(gates, reverse(gates))

# for t in 0.0:tau:ttotal
#     flux = mean(calculate_spinFlux(Ivec, rho, sites, [1:L;]))
#     println("$t $flux")
#     # push!(fluxVec, real.(flux))

#     t≈ttotal && break
    
#     # global rho = apply(gates, rho; cutoff=cutoff)

#     global rho = product(gates, rho)
#     global rho = truncate(rho, cutoff=cutoff)

#     # normalize!(rho)
#     # trRho = rho * delta(i,l)
#     # println(tr(ITensor(rho)))
# end

currents = calculate_spinFlux(Ivec, rho, sites, [1:L;])
println("Analytical value of the current: ", -fluxXXX(L, gamma))
println("Numerical approximation of the current: ", mean(currents))
# println("Experimentation: ", mean(fluxVec))

analyticalMag = magnonDensityXXX(L, gamma)
magDensity = calculate_magnetization(Ivec, rho, sites)

hline!([-fluxXXX(L, gamma)], linewidth=2, color=:black, linestyle=:solid, label="Analytical value")
plot!(fluxVec)
# hline!([-fluxXXX(L, gamma)], linewidth=2, color=:black, linestyle=:solid, label="Analytical value")

# println("Analytical Magnetization profile: ", analyticalMag)
# println("Numerical magProf Approx (real) : ", real.(magDensity))
# println("Numerical magProf Approx (cmpx) : ", imag.(magDensity))
# plot([1:L;], real.(analyticalMag))
# Using the imaginary part of the magnon density as  uncertainty (clearly an overestimatation)
# plot!([1:L;], real.(magDensity),yerr=imag.(magDensity), linestyle=:dash)
# annotate!(1.5, -0.15, text(string("Analytical value of the current: ", -fluxXXX(L, gamma)), 10))
ylims!(-0.6, 0)
gui()

# savefig("timeEvolutionFromDMRGansatz.pdf")
savefig("timeEvolutionFromDMRGGS.pdf")


# Keeps plots open
read(stdin, Char)
