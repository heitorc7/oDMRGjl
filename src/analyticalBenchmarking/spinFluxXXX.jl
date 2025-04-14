using LinearAlgebra
using SparseArrays

# As derived by G.T. Landi and D. Karevski in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.174422
# Usage: flux(N-> system size, gamma-> coupling parameter with the baths)
function fluxXXX(N::Int, γγ::Float64; chop_tol=1e-10)
    γ = γγ / 8.0
    h = 0.0
    p = im / (2.0 * (γ - im * h))
    q = -im / (2.0 * (γ + im * h))
    
    # Exact Mathematica-style matrix construction
    B0 = spzeros(ComplexF64, N+1, N+1)
    for k in 1:N+1
        for l in 1:N+1
            if k == l
                B0[k,l] = 2*(p-(k-1))*(q-(k-1))
            elseif k == l-1
                B0[k,l] = (l-1)^2
            elseif k == l+1
                B0[k,l] = (2p-(l-1))*(2q-(l-1))
            end
        end
    end
    
    v = zeros(ComplexF64, N+1)
    v[1] = 1.0  # NEye[][[1]] equivalent
    
    # Stabilized matrix power-vector product
    tmp = copy(v)
    for _ in 1:N-1
        tmp = B0 * tmp
        tmp ./= norm(tmp)  # Normalize at each step
    end
    
    # Exact Mathematica computation order
    num = 4im * (q - p) * dot(v, tmp)
    den = dot(v, B0 * tmp)
    result = num / (2 * den)
    
    # Mathematica's Chop equivalent
    abs(imag(result)) < chop_tol || @warn "Significant imaginary part"
    real_part = real(result)
    return abs(real_part) < chop_tol ? 0.0 : real_part
end