using LinearAlgebra
using SparseArrays

# As derived by G.T. Landi and D. Karevski in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.174422
# Usage: magnonDensity(N-> system size, gamma-> coupling parameter with the baths)
function magnonDensityXXX(N::Int, γγ::Float64, δi::Int=1; chop_tol=1e-10)
    h = 0.0
    γ = γγ / 8.0
    p = im / (2 * (γ - im * h))
    
    # Construct B0 matrix
    B0 = spzeros(ComplexF64, N+1, N+1)
    for k in 1:N+1
        for l in 1:N+1
            if k == l
                B0[k,l] = 2*(p - (k-1)) * conj(p - (k-1))
            elseif k == l-1
                B0[k,l] = (l-1)^2
            elseif k == l+1
                B0[k,l] = (2p - (l-1)) * conj(2p - (l-1))
            end
        end
    end
    
    # Construct Bz matrix (anti-Hermitian part)
    Bz = spzeros(ComplexF64, N+1, N+1)
    for k in 1:N
        Bz[k,k+1] = k^2
        Bz[k+1,k] = -abs(2p - k)^2  # Note sign flip from B0
    end
    
    v = zeros(ComplexF64, N+1)
    v[1] = 1.0  # NEye[][[1]] equivalent
    
    results = Vector{Float64}(undef, length(1:δi:N))
    
    for (idx, i) in enumerate(1:δi:N)
        # Stabilized matrix powers with normalization
        Iv = v
        if i > 1
            Iv = v
            for _ in 1:i-1
                Iv = B0' * Iv
                Iv ./= norm(Iv)
            end
        end
        
        Rv = v
        if i < N
            Rv = v
            for _ in 1:N-i
                Rv = B0 * Rv
                Rv ./= norm(Rv)
            end
        end
        
        # Compute the ratio with stabilization
        numerator = dot(Iv, Bz * Rv)
        denominator = dot(Iv, B0 * Rv)
        ratio = numerator / denominator
        
        # Mathematica's exact formula with chop
        val = -2 * 0.5 * (1 + ratio) + 1
        results[idx] = abs(imag(val)) < chop_tol ? real(val) : error("Significant imaginary part detected")
    end
    
    return results
end