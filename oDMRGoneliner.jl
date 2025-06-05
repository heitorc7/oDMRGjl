using ITensors
using ITensorMPS
using Plots
using Statistics
include("src/sites/twoSpinHalf.jl")
include("src/objects/liouvillian.jl")
include("src/objects/Ivec.jl")
include("src/observables/spinFlux.jl")
include("src/observables/magnonDensity.jl")
include("src/funcs/warmUp.jl")
include("src/funcs/fileHandling.jl")
include("src/funcs/oDMRGobserver.jl")
# include("src/analyticalBenchmarking/spinFluxXXX.jl")
# include("src/analyticalBenchmarking/magnonDensityXXX.jl")

function main(args)
    # Parse command line arguments
    length(args) == 4 || error("Usage: julia oneLiner.jl N gamma J Delta")

    # Physical parameters: system size and coupling factor gamma
    L = parse(Int, args[1])
    gamma = parse(Float64, args[2])
    J = parse(Float64, args[3])
    Delta = parse(Float64, args[4])

    # Physical Parameters: Jx, Jy, Jz, dissipators location, dissipators temperature, coupling vector, magnetic field
    params = (
        Jx = J, Jy = J, Jz = Delta,
        dissipatorsVec = [1, L],
        dissipatorsTempValues = [1.0, 0.0],
        gammaVec = [gamma, gamma],
        hVec = zeros(L)
    )
    # Create sites using customized TwoSpinHalf sites
    sites = siteinds("TwoSpinHalf", L)

    save = true
    global iter = 1

    # Create output directory structure
    output_dir = "simulationData/N_$(L)/"
    mkpath(output_dir)  # Creates directory if not exists
    # Generate output filename (.dat recommended for numerical data)
    output_file = joinpath(output_dir, "XXZ_N-$(L)_gamma-$(gamma)_Jxy-$(J)_Delta-$(Delta)_iter-$(iter).txt")
    while isfile(output_file)
        println("Found a file with desired parameter set, incresing iter number to not override it. Current iter number: $(iter)")
        global iter += 1
        output_file = joinpath(output_dir, "XXZ_N-$(L)_gamma-$(gamma)_Jxy-$(J)_Delta-$(Delta)_iter-$(iter).txt")
    end

    open(output_file, "w") do file
        redirect_stdout(file) do

            # LdL MPO constructor
            global LdL = LdLXYZConstruct(sites, params...; maxdim=500, cutoff=1e-8)
            global rho = make_ivec(L, sites)

            cutoff = 1E-10

            warmUp(rho, LdL, 0.0000001)
            global Ivec = make_ivec(L, sites)

            # Now initiating iterative maxDim increase oDMRG
            global maxdim2 = 2;
            maxDimInc = 2;
            sweepBDChangeThresholdValue = 0.001 # good values: 0.0005 ~ 0.001
            finalSweep = 1000; # usually 1000~2000 works for N = 10. Experiment with more for larger sizes

            global energyFin = floatmax(Float64)
            global energyThreshold = 1E-7
            obs = oDMRGobserver(1E-7)

            global sweep = 1
            while abs(energyFin) > energyThreshold
                global sweep += 1
                # for sweep in 1:finalSweep
                energyIni = energyFin
                global energyFin, rho = dmrg(LdL,rho; nsweeps=1, maxdim=maxdim2, cutoff=cutoff)

                # Now we calculate the flux
                # currents = calculate_spinFlux(Ivec, rho, sites, [1:L;])
                if (abs(energyFin) < abs(energyIni)) && abs(abs(energyFin) - abs(energyIni))/energyIni < sweepBDChangeThresholdValue
                    global maxdim2 += maxDimInc
                end

                if sweep%10 == 0
                    # Calculating observables: currents and magDensity for every site
                    global currents = calculate_spinFlux(Ivec, rho, sites, [1:L;])
                    global magDensity = calculate_magnetization(Ivec, rho, sites)

                    println(currents)
                    println(magDensity)
                    if (sweep%500 == 0) & save
                        saveMPS(rho, "myMPS")
                    end
                end
            end
            println("Energy under the desired threshold of $energyThreshold. Final energy = $energyFin")

            println("Numerical approximation of the current: ", currents)
            println("Numerical mean of the currents: ", mean(currents))

            println("Numerical magProf Approx (real) : ", real.(magDensity))
            println("Numerical magProf Approx (cmpx) : ", imag.(magDensity))

            # currents = calculate_spinFlux(Ivec, rho, sites, [1:L;])
            # println("Analytical value of the current: ", -fluxXXX(L, gamma))
            # println("Numerical approximation of the current: ", mean(currents))

            # analyticalMag = magnonDensityXXX(L, gamma)
            # magDensity = calculate_magnetization(Ivec, rho, sites)

            # println("Analytical Magnetization profile: ", analyticalMag)
            # println("Numerical magProf Approx (real) : ", real.(magDensity))
            # println("Numerical magProf Approx (cmpx) : ", imag.(magDensity))
            # plot([1:L;], real.(analyticalMag))
            # # Using the imaginary part of the magnon density as  uncertainty (clearly an overestimatation)
            # plot!([1:L;], real.(magDensity),yerr=imag.(magDensity), linestyle=:dash)
            # gui()

            # # Keeps plots open
            # read(stdin, Char)

            if save
                saveMPS(rho, "myMPS")
            end
        end
    end
end

# Execute main
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end