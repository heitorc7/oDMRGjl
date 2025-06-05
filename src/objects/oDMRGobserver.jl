using ITensors, ITensorMPS

mutable struct oDMRGobserver <: AbstractObserver
   energy_tol::Float64
   last_energy::Float64

   oDMRGobserver(energy_tol=1E-7) = new(energy_tol,1000.0)
end

function ITensorMPS.checkdone!(o::oDMRGobserver;kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    if abs(energy-o.last_energy)/abs(energy) < o.energy_tol
      println("Stopping DMRG after sweep $sw")
      return true
    end
    # Otherwise, update last_energy and keep going
    o.last_energy = energy
    return false
  end

  function ITensorMPS.measure!(o::oDMRGobserver; kwargs...)
    energy = kwargs[:energy]
    sweep = kwargs[:sweep]
    bond = kwargs[:bond]
    outputlevel = kwargs[:outputlevel]
  
    # if outputlevel > 0
    #   println("Sweep $sweep at bond $bond, the energy is $energy")
    # end
  end